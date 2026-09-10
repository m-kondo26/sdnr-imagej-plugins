import ij.IJ;
import ij.ImagePlus;
import ij.WindowManager;
import ij.gui.GenericDialog;
import ij.gui.NonBlockingGenericDialog;
import ij.gui.Overlay;
import ij.gui.Plot;
import ij.gui.PlotWindow;
import ij.gui.Roi;
import ij.io.FileInfo;
import ij.io.OpenDialog;
import ij.io.RoiDecoder;
import ij.io.SaveDialog;
import ij.measure.Calibration;
import ij.measure.ResultsTable;
import ij.plugin.PlugIn;
import ij.process.FHT;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;

import java.awt.Color;
import java.awt.Rectangle;
import java.awt.Window;
import java.io.ByteArrayOutputStream;
import java.io.FileInputStream;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.List;
import java.util.Locale;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.zip.ZipEntry;
import java.util.zip.ZipInputStream;
import java.util.zip.ZipOutputStream;

/**
 * Batch NPS measurement driven by multiple ImageJ ROI zip entries.
 */
public class NPS_2D_FFT_ROI_Batch implements PlugIn {
    private static final String VERSION = "2026-09-09-nps-input-validation";
    private static final int FIRST_OUTPUT_BIN = 1;
    private static final String TREND_REMOVAL_QUADRATIC = "2D second-order polynomial (Samei)";
    private static final String TREND_REMOVAL_MEAN = "Subtract ROI mean only";
    private static final String TREND_REMOVAL_NONE = "None";
    private static final String[] TREND_REMOVAL_CHOICES = {
        TREND_REMOVAL_QUADRATIC,
        TREND_REMOVAL_MEAN,
        TREND_REMOVAL_NONE
    };
    private static final Pattern NUMBER_PATTERN = Pattern.compile("[-+]?\\d*\\.?\\d+(?:[Ee][-+]?\\d+)?");

    @Override
    public void run(String arg) {
        IJ.log("NPS 2D FFT ROI Batch " + VERSION + " loaded as class " + getClass().getName());
        ImagePlus imp = IJ.getImage();
        if (imp == null) {
            IJ.error("No image", "Open an image stack first.");
            return;
        }
        if (imp.getType() == ImagePlus.COLOR_RGB) {
            IJ.error("Unsupported image", "NPS should be calculated from a grayscale CT image.");
            return;
        }
        if (imp.getStackSize() < 1) {
            IJ.error("No stack", "The current image has no slices.");
            return;
        }

        Calibration cal = imp.getCalibration();
        Geometry geometry;
        try {
            geometry = inferGeometry(imp, cal, imp.getWidth(), imp.getHeight());
        } catch (IllegalArgumentException e) {
            IJ.error("Invalid NPS geometry", e.getMessage());
            return;
        }
        BatchParameters params = readBatchParameters(imp, geometry);
        if (params == null) {
            return;
        }

        String roiZipPath = chooseRoiZip();
        if (roiZipPath == null) {
            return;
        }

        List<RoiEntry> roiEntries;
        try {
            roiEntries = readRoiZip(roiZipPath);
        } catch (IOException e) {
            IJ.error("ROI zip read failed", e.getMessage());
            return;
        }
        if (roiEntries.isEmpty()) {
            IJ.error("No ROI", "The selected zip file does not contain readable ImageJ .roi entries.");
            return;
        }
        assignUniqueDisplayNames(roiEntries);

        List<RoiEntry> validRois = validateRois(roiEntries, imp.getWidth(), imp.getHeight(), params.matrix);
        if (validRois == null) {
            clearStatusAndProgress();
            return;
        }
        if (!allSameRoiSize(validRois)) {
            IJ.error(
                "Mixed ROI sizes",
                "All valid ROIs must have the same square size to average NPS directly. Use same-size ROI entries in the ROI zip."
            );
            clearStatusAndProgress();
            return;
        }

        showMeasuredRoiPreview(imp, validRois);
        NpsResult result;
        try {
            result = computeRoiAveragedNps(imp, cal, params, validRois);
        } catch (IllegalArgumentException e) {
            IJ.error("NPS calculation stopped", e.getMessage() + "\nNo workbook was exported.");
            clearStatusAndProgress();
            return;
        }
        showSummary(imp, params, validRois, result);
        showNpsPlot(result);

        boolean saved = exportWorkbook(params, result);
        IJ.log(String.format(
            Locale.US,
            "NPS ROI Batch finished: image=%s rois=%d slices=%d export=%s",
            imp.getTitle(),
            validRois.size(),
            params.sliceNumbers.length,
            saved ? "saved" : "not saved"
        ));
        clearStatusAndProgress();
    }

    private BatchParameters readBatchParameters(ImagePlus imp, Geometry geometry) {
        int currentSlice = Math.max(1, Math.min(imp.getCurrentSlice(), imp.getStackSize()));
        GenericDialog gd = new NonBlockingGenericDialog("NPS ROI Batch Parameters");
        gd.addMessage(String.format(
            Locale.US,
            "Geometry: FOV %.9g mm / matrix %d = %.9g mm/pixel%nSource: %s",
            geometry.fovMm,
            geometry.matrix,
            geometry.pixelSizeMm,
            geometry.source
        ));
        if (!geometry.isDicomDerived) {
            gd.addMessage("DICOM FOV/pixel-spacing tags were not found. Confirm or edit FOV and matrix below.");
        }
        gd.addMessage("After this dialog, select an ImageJ ROI zip. All valid rectangular square power-of-two ROIs will be measured.");
        gd.addNumericField("FOV (mm)", geometry.fovMm, 6);
        gd.addNumericField("Matrix size (pixels)", geometry.matrix, 0);
        gd.addNumericField("Start slice", currentSlice, 0);
        gd.addNumericField("End slice", imp.getStackSize(), 0);
        gd.addNumericField("Slice step", 1, 0);
        gd.addChoice("Trend removal before FFT", TREND_REMOVAL_CHOICES, TREND_REMOVAL_QUADRATIC);
        gd.addCheckbox("Apply 2D Hann window", false);
        gd.showDialog();
        if (gd.wasCanceled()) {
            return null;
        }

        double enteredFov = gd.getNextNumber();
        double enteredMatrix = gd.getNextNumber();
        double enteredStart = gd.getNextNumber();
        double enteredEnd = gd.getNextNumber();
        double enteredStep = gd.getNextNumber();
        String trendRemoval = gd.getNextChoice();
        boolean useHann = gd.getNextBoolean();
        if (gd.invalidNumber()
            || enteredFov <= 0 || enteredMatrix <= 0
            || enteredStart < 1 || enteredEnd < 1 || enteredStep < 1
            || !Double.isFinite(enteredFov) || !Double.isFinite(enteredMatrix)
            || !Double.isFinite(enteredStart) || !Double.isFinite(enteredEnd) || !Double.isFinite(enteredStep)) {
            IJ.error("Invalid input", "FOV, matrix, and slice range must be valid positive numbers.");
            return null;
        }

        int matrix = Math.max(1, (int) Math.round(enteredMatrix));
        int startSlice = clampSlice((int) Math.round(enteredStart), imp.getStackSize());
        int endSlice = clampSlice((int) Math.round(enteredEnd), imp.getStackSize());
        if (endSlice < startSlice) {
            IJ.error("Invalid slice range", "End slice must be greater than or equal to start slice.");
            return null;
        }
        int step = Math.max(1, (int) Math.round(enteredStep));
        double pixelSizeMm = enteredFov / matrix;
        boolean unchangedGeometry = nearlyEqual(enteredFov, geometry.fovMm) && matrix == geometry.matrix;
        String source = unchangedGeometry
            ? geometry.source
            : String.format(Locale.US, "Manual dialog input: FOV=%.9g mm, matrix=%d", enteredFov, matrix);
        return new BatchParameters(
            pixelSizeMm,
            enteredFov,
            matrix,
            geometry.isDicomDerived && unchangedGeometry,
            source,
            sliceNumbers(startSlice, endSlice, step),
            startSlice,
            endSlice,
            step,
            trendRemoval,
            useHann
        );
    }

    private String chooseRoiZip() {
        OpenDialog od = new OpenDialog("Open ImageJ ROI zip", null);
        String fileName = od.getFileName();
        String directory = od.getDirectory();
        if (fileName == null || directory == null) {
            IJ.log("NPS ROI Batch cancelled: ROI zip was not selected.");
            return null;
        }
        return Paths.get(directory, fileName).toString();
    }

    private List<RoiEntry> readRoiZip(String path) throws IOException {
        List<RoiEntry> entries = new ArrayList<>();
        byte[] buffer = new byte[8192];
        ZipInputStream zis = new ZipInputStream(new FileInputStream(path));
        try {
            ZipEntry entry;
            while ((entry = zis.getNextEntry()) != null) {
                if (entry.isDirectory() || !entry.getName().toLowerCase(Locale.US).endsWith(".roi")) {
                    continue;
                }
                ByteArrayOutputStream bytes = new ByteArrayOutputStream();
                int read;
                while ((read = zis.read(buffer)) >= 0) {
                    bytes.write(buffer, 0, read);
                }
                Roi roi = new RoiDecoder(bytes.toByteArray(), entry.getName()).getRoi();
                if (roi == null) {
                    IJ.log("NPS ROI Batch: skipped unreadable ROI entry " + entry.getName());
                    continue;
                }
                String name = roi.getName();
                if (name == null || name.trim().isEmpty()) {
                    name = baseName(entry.getName());
                    roi.setName(name);
                }
                entries.add(new RoiEntry(roi, name, entry.getName()));
            }
        } finally {
            zis.close();
        }
        return entries;
    }

    private void assignUniqueDisplayNames(List<RoiEntry> entries) {
        for (int i = 0; i < entries.size(); i++) {
            RoiEntry current = entries.get(i);
            boolean duplicate = false;
            for (int j = 0; j < entries.size(); j++) {
                if (i != j && current.name.equals(entries.get(j).name)) {
                    duplicate = true;
                    break;
                }
            }
            current.displayName = duplicate
                ? current.name + " [" + current.entryName + "]"
                : current.name;
        }
    }

    private List<RoiEntry> validateRois(List<RoiEntry> entries, int width, int height, int matrixLimit) {
        List<String> invalid = new ArrayList<>();
        for (RoiEntry entry : entries) {
            String reason = invalidRoiReason(entry.roi, width, height, matrixLimit);
            if (reason == null) {
                continue;
            }
            String message = entry.displayName + ": " + reason;
            invalid.add(message);
            IJ.log("NPS ROI Batch: invalid ROI " + message);
        }
        if (!invalid.isEmpty()) {
            StringBuilder message = new StringBuilder();
            message.append("All ROI entries must be rectangular square ROIs whose side length is a power of two, inside the image, and not larger than matrix size.\n\n");
            int shown = Math.min(invalid.size(), 8);
            for (int i = 0; i < shown; i++) {
                message.append("- ").append(invalid.get(i)).append('\n');
            }
            if (invalid.size() > shown) {
                message.append("- ... ").append(invalid.size() - shown).append(" more invalid ROI entries");
            } else {
                int last = message.length() - 1;
                if (last >= 0 && message.charAt(last) == '\n') {
                    message.deleteCharAt(last);
                }
            }
            IJ.error("Invalid ROI", message.toString());
            return null;
        }
        IJ.log(String.format(
            Locale.US,
            "NPS ROI Batch: accepted %d/%d ROI entries from zip.",
            entries.size(),
            entries.size()
        ));
        return new ArrayList<>(entries);
    }

    private String invalidRoiReason(Roi roi, int width, int height, int matrixLimit) {
        if (roi == null) {
            return "ROI is null";
        }
        if (roi.getType() != Roi.RECTANGLE || roi.getCornerDiameter() != 0) {
            return "ROI must be a rectangle without rounded corners";
        }
        Rectangle bounds = roi.getBounds();
        if (bounds.width < 8 || bounds.height < 8) {
            return "ROI must be at least 8x8 pixels";
        }
        if (bounds.width != bounds.height) {
            return "ROI is not square";
        }
        if (!isPowerOfTwo(bounds.width)) {
            return "ROI side length is not a power of two";
        }
        if (bounds.width > matrixLimit) {
            return "ROI side length is larger than matrix size";
        }
        if (bounds.x < 0 || bounds.y < 0 || bounds.x + bounds.width > width || bounds.y + bounds.height > height) {
            return "ROI is outside the image";
        }
        return null;
    }

    private boolean allSameRoiSize(List<RoiEntry> rois) {
        if (rois.isEmpty()) {
            return true;
        }
        int expected = rois.get(0).roi.getBounds().width;
        for (RoiEntry roi : rois) {
            int size = roi.roi.getBounds().width;
            if (size != expected) {
                IJ.log(String.format(
                    Locale.US,
                    "NPS ROI Batch: mixed valid ROI sizes are not averaged directly. First ROI size=%d px, ROI %s size=%d px.",
                    expected,
                    roi.displayName,
                    size
                ));
                return false;
            }
        }
        return true;
    }

    private NpsResult computeRoiAveragedNps(ImagePlus imp, Calibration cal, BatchParameters params, List<RoiEntry> rois) {
        if (rois == null || rois.isEmpty()) {
            throw new IllegalArgumentException("At least one valid NPS ROI is required.");
        }
        for (RoiEntry roi : rois) {
            String reason = invalidRoiReason(roi.roi, imp.getWidth(), imp.getHeight(), params.matrix);
            if (reason != null) {
                throw new IllegalArgumentException("ROI " + roi.displayName + ": " + reason);
            }
        }
        if (!allSameRoiSize(rois)) {
            throw new IllegalArgumentException("All averaged NPS ROIs must have the same size.");
        }
        NpsResult first = null;
        double[][] roiAveragedPerSlice = null;

        for (int i = 0; i < rois.size(); i++) {
            RoiEntry roi = rois.get(i);
            IJ.showStatus("NPS ROI Batch: processing " + roi.displayName + " (" + (i + 1) + "/" + rois.size() + ")");
            IJ.showProgress(i, rois.size());
            NpsResult roiResult;
            try {
                roiResult = computeNps(imp, cal, params.forRoi(roi));
            } catch (IllegalArgumentException e) {
                throw new IllegalArgumentException("ROI " + roi.displayName + ": " + e.getMessage(), e);
            }
            if (first == null) {
                first = roiResult;
                roiAveragedPerSlice = new double[roiResult.perSlice.length][roiResult.frequencies.length];
            }
            for (int s = 0; s < roiResult.perSlice.length; s++) {
                for (int r = 0; r < roiResult.frequencies.length; r++) {
                    roiAveragedPerSlice[s][r] += roiResult.perSlice[s][r];
                }
            }
        }

        for (int s = 0; s < roiAveragedPerSlice.length; s++) {
            for (int r = 0; r < roiAveragedPerSlice[s].length; r++) {
                roiAveragedPerSlice[s][r] /= rois.size();
                if (!Double.isFinite(roiAveragedPerSlice[s][r])) {
                    throw new IllegalArgumentException("Nonfinite ROI-averaged NPS at frequency bin " + r + ".");
                }
            }
        }

        double[] sliceAverage = new double[first.frequencies.length];
        for (int r = 0; r < sliceAverage.length; r++) {
            for (int s = 0; s < roiAveragedPerSlice.length; s++) {
                sliceAverage[r] += roiAveragedPerSlice[s][r];
            }
            sliceAverage[r] /= roiAveragedPerSlice.length;
            if (!Double.isFinite(sliceAverage[r])) {
                throw new IllegalArgumentException("Nonfinite average NPS at frequency bin " + r + ".");
            }
        }

        return new NpsResult(
            first.frequencies,
            roiAveragedPerSlice,
            sliceAverage,
            first.sliceNumbers,
            first.frequencyStep,
            first.nyquist,
            first.parameters,
            rois.size()
        );
    }

    private void validateNpsInput(ImagePlus imp, NpsParameters params) {
        if (imp == null || params == null) {
            throw new IllegalArgumentException("An image and NPS parameters are required.");
        }
        double pixel = params.pixelSizeMm;
        if (!Double.isFinite(pixel) || pixel <= 0 || !Double.isFinite(pixel * pixel) || pixel * pixel <= 0) {
            throw new IllegalArgumentException("Pixel size must be finite, positive, and representable in mm.");
        }
        Rectangle roi = params.roiBounds;
        if (roi == null || roi.width < 8 || roi.width != roi.height || !isPowerOfTwo(roi.width)
            || roi.width > params.matrix || roi.x < 0 || roi.y < 0
            || (long) roi.x + roi.width > imp.getWidth() || (long) roi.y + roi.height > imp.getHeight()) {
            throw new IllegalArgumentException("NPS ROI must be a square of at least 8 pixels, a power of two, within the image and matrix.");
        }
        double step = 1.0 / (roi.width * pixel);
        if (!Double.isFinite(step) || step <= 0) {
            throw new IllegalArgumentException("The ROI size and pixel size produce an invalid frequency step.");
        }
        if (params.sliceNumbers == null || params.sliceNumbers.length == 0) {
            throw new IllegalArgumentException("At least one slice is required.");
        }
        for (int slice : params.sliceNumbers) {
            if (slice < 1 || slice > imp.getStackSize()) {
                throw new IllegalArgumentException("Slice " + slice + " is outside the image stack.");
            }
        }
    }

    private NpsResult computeNps(ImagePlus imp, Calibration cal, NpsParameters params) {
        validateNpsInput(imp, params);
        int roiSize = params.roiBounds.width;
        int maxRadius = roiSize / 2;
        double pixelAreaMm2 = params.pixelSizeMm * params.pixelSizeMm;
        double frequencyStep = 1.0 / (roiSize * params.pixelSizeMm);
        double[] frequencies = new double[maxRadius + 1];
        for (int r = 0; r <= maxRadius; r++) {
            frequencies[r] = r * frequencyStep;
        }

        double[] window = buildHannWindow(roiSize, params.useHannWindow);
        double meanWindowSquared = mean2DWindowSquared(window);
        double[][] perSlice = new double[params.sliceNumbers.length][maxRadius + 1];
        double[] average = new double[maxRadius + 1];

        for (int i = 0; i < params.sliceNumbers.length; i++) {
            int slice = params.sliceNumbers[i];
            FloatProcessor fp;
            try {
                fp = calibratedRoiToFloat(imp.getStack().getProcessor(slice), cal, params.roiBounds);
            } catch (IllegalArgumentException e) {
                throw new IllegalArgumentException("Slice " + slice + ": " + e.getMessage(), e);
            }
            applyTrendRemoval(fp, params.trendRemoval);
            applySeparableWindow(fp, window);

            FHT fht = new FHT(fp);
            fht.transform();
            FloatProcessor power = fht.getRawPowerSpectrum();
            fht.swapQuadrants(power);

            double[] radial = radialAverage(power, maxRadius);
            for (int r = 0; r <= maxRadius; r++) {
                double nps = radial[r] * pixelAreaMm2 / (roiSize * roiSize) / meanWindowSquared;
                if (!Double.isFinite(nps) || nps < 0) {
                    throw new IllegalArgumentException("Slice " + slice + ": nonfinite or negative NPS at frequency bin " + r + ".");
                }
                perSlice[i][r] = nps;
                average[r] += nps;
            }
        }

        for (int r = 0; r <= maxRadius; r++) {
            average[r] /= params.sliceNumbers.length;
            if (!Double.isFinite(average[r])) {
                throw new IllegalArgumentException("Nonfinite average NPS at frequency bin " + r + ".");
            }
        }
        return new NpsResult(
            frequencies,
            perSlice,
            average,
            params.sliceNumbers,
            frequencyStep,
            1.0 / (2.0 * params.pixelSizeMm),
            params,
            1
        );
    }

    private FloatProcessor calibratedRoiToFloat(ImageProcessor ip, Calibration cal, Rectangle roi) {
        FloatProcessor fp = new FloatProcessor(roi.width, roi.height);
        for (int y = 0; y < roi.height; y++) {
            for (int x = 0; x < roi.width; x++) {
                double raw = ip.getf(roi.x + x, roi.y + y);
                double value = cal != null ? cal.getCValue(raw) : raw;
                if (!Double.isFinite(raw) || !Double.isFinite(value) || !Float.isFinite((float) value)) {
                    throw new IllegalArgumentException("Nonfinite or out-of-range intensity at image pixel ("
                        + (roi.x + x) + ", " + (roi.y + y) + ").");
                }
                fp.setf(x, y, (float) value);
            }
        }
        return fp;
    }

    private void applyTrendRemoval(FloatProcessor fp, String trendRemoval) {
        if (TREND_REMOVAL_QUADRATIC.equals(trendRemoval)) {
            subtractSecondOrderPolynomialTrend(fp);
        } else if (TREND_REMOVAL_MEAN.equals(trendRemoval)) {
            subtractMean(fp);
        }
    }

    private void subtractMean(FloatProcessor fp) {
        float[] pixels = (float[]) fp.getPixels();
        double sum = 0.0;
        for (float pixel : pixels) {
            sum += pixel;
        }
        double mean = sum / pixels.length;
        for (int i = 0; i < pixels.length; i++) {
            pixels[i] = (float) (pixels[i] - mean);
        }
    }

    private void subtractSecondOrderPolynomialTrend(FloatProcessor fp) {
        double[] coeff = fitSecondOrderPolynomial(fp);
        if (coeff == null) {
            IJ.log("NPS ROI Batch: 2D second-order trend fit failed; falling back to ROI mean subtraction.");
            subtractMean(fp);
            return;
        }

        int width = fp.getWidth();
        int height = fp.getHeight();
        float[] pixels = (float[]) fp.getPixels();
        double xCenter = (width - 1) / 2.0;
        double yCenter = (height - 1) / 2.0;
        double xScale = width > 1 ? 2.0 / (width - 1) : 1.0;
        double yScale = height > 1 ? 2.0 / (height - 1) : 1.0;

        for (int y = 0; y < height; y++) {
            double yn = (y - yCenter) * yScale;
            double y2 = yn * yn;
            for (int x = 0; x < width; x++) {
                double xn = (x - xCenter) * xScale;
                double fit = coeff[0]
                    + coeff[1] * xn
                    + coeff[2] * yn
                    + coeff[3] * xn * xn
                    + coeff[4] * xn * yn
                    + coeff[5] * y2;
                int index = y * width + x;
                pixels[index] = (float) (pixels[index] - fit);
            }
        }
    }

    private double[] fitSecondOrderPolynomial(FloatProcessor fp) {
        int width = fp.getWidth();
        int height = fp.getHeight();
        float[] pixels = (float[]) fp.getPixels();
        double[][] normal = new double[6][6];
        double[] rhs = new double[6];
        double[] basis = new double[6];
        double xCenter = (width - 1) / 2.0;
        double yCenter = (height - 1) / 2.0;
        double xScale = width > 1 ? 2.0 / (width - 1) : 1.0;
        double yScale = height > 1 ? 2.0 / (height - 1) : 1.0;

        for (int y = 0; y < height; y++) {
            double yn = (y - yCenter) * yScale;
            for (int x = 0; x < width; x++) {
                double xn = (x - xCenter) * xScale;
                basis[0] = 1.0;
                basis[1] = xn;
                basis[2] = yn;
                basis[3] = xn * xn;
                basis[4] = xn * yn;
                basis[5] = yn * yn;
                double value = pixels[y * width + x];
                for (int i = 0; i < basis.length; i++) {
                    rhs[i] += basis[i] * value;
                    for (int j = i; j < basis.length; j++) {
                        normal[i][j] += basis[i] * basis[j];
                    }
                }
            }
        }

        for (int i = 0; i < normal.length; i++) {
            for (int j = 0; j < i; j++) {
                normal[i][j] = normal[j][i];
            }
        }
        return solveLinearSystem(normal, rhs);
    }

    private double[] solveLinearSystem(double[][] matrix, double[] rhs) {
        int n = rhs.length;
        double[][] a = new double[n][n + 1];
        for (int row = 0; row < n; row++) {
            System.arraycopy(matrix[row], 0, a[row], 0, n);
            a[row][n] = rhs[row];
        }

        for (int col = 0; col < n; col++) {
            int pivot = col;
            double max = Math.abs(a[col][col]);
            for (int row = col + 1; row < n; row++) {
                double value = Math.abs(a[row][col]);
                if (value > max) {
                    max = value;
                    pivot = row;
                }
            }
            if (max < 1e-12) {
                return null;
            }
            if (pivot != col) {
                double[] tmp = a[col];
                a[col] = a[pivot];
                a[pivot] = tmp;
            }

            for (int row = col + 1; row < n; row++) {
                double factor = a[row][col] / a[col][col];
                a[row][col] = 0.0;
                for (int c = col + 1; c <= n; c++) {
                    a[row][c] -= factor * a[col][c];
                }
            }
        }

        double[] solution = new double[n];
        for (int row = n - 1; row >= 0; row--) {
            double sum = a[row][n];
            for (int col = row + 1; col < n; col++) {
                sum -= a[row][col] * solution[col];
            }
            solution[row] = sum / a[row][row];
        }
        return solution;
    }

    private double[] buildHannWindow(int size, boolean useHann) {
        double[] window = new double[size];
        if (!useHann) {
            for (int i = 0; i < size; i++) {
                window[i] = 1.0;
            }
            return window;
        }
        for (int i = 0; i < size; i++) {
            window[i] = 0.5 * (1.0 - Math.cos(2.0 * Math.PI * i / (size - 1)));
        }
        return window;
    }

    private double mean2DWindowSquared(double[] window) {
        double sumW2 = 0.0;
        for (double w : window) {
            sumW2 += w * w;
        }
        return (sumW2 * sumW2) / (window.length * window.length);
    }

    private void applySeparableWindow(FloatProcessor fp, double[] window) {
        int size = window.length;
        for (int y = 0; y < size; y++) {
            double wy = window[y];
            for (int x = 0; x < size; x++) {
                fp.setf(x, y, (float) (fp.getf(x, y) * window[x] * wy));
            }
        }
    }

    private double[] radialAverage(FloatProcessor power, int maxRadius) {
        int size = power.getWidth();
        double[] sums = new double[maxRadius + 1];
        int[] counts = new int[maxRadius + 1];
        int center = size / 2;
        float[] values = (float[]) power.getPixels();
        for (int y = 0; y < size; y++) {
            int dy = y - center;
            for (int x = 0; x < size; x++) {
                int dx = x - center;
                int r = (int) Math.floor(Math.hypot(dx, dy));
                if (r <= maxRadius) {
                    sums[r] += values[y * size + x];
                    counts[r]++;
                }
            }
        }

        double[] average = new double[maxRadius + 1];
        for (int r = 0; r <= maxRadius; r++) {
            average[r] = counts[r] > 0 ? sums[r] / counts[r] : Double.NaN;
        }
        return average;
    }

    private void showSummary(ImagePlus imp, BatchParameters params, List<RoiEntry> rois, NpsResult result) {
        ResultsTable rt = ResultsTable.getResultsTable();
        if (rt == null) {
            rt = new ResultsTable();
        }
        rt.incrementCounter();
        rt.addValue("ROI_Count", rois.size());
        rt.addValue("Slices", result.perSlice.length);
        rt.addValue("PixelSize_mm", params.pixelSizeMm);
        rt.addValue("FOV_mm", params.fovMm);
        rt.addValue("Matrix", params.matrix);
        rt.addValue("Geometry_DICOM_Derived", params.geometryDicomDerived ? 1 : 0);
        rt.addValue("ROI_Size_px", rois.get(0).roi.getBounds().width);
        rt.addValue("Quadratic_Trend_Removed", TREND_REMOVAL_QUADRATIC.equals(params.trendRemoval) ? 1 : 0);
        rt.addValue("Mean_Subtracted", TREND_REMOVAL_NONE.equals(params.trendRemoval) ? 0 : 1);
        rt.addValue("Hann_Window_Applied", params.useHannWindow ? 1 : 0);
        rt.addValue("FFT_Length_mm", rois.get(0).roi.getBounds().width * params.pixelSizeMm);
        rt.addValue("FrequencyStep_cycles_per_mm", result.frequencyStep);
        rt.addValue("Nyquist_cycles_per_mm", result.nyquist);
        rt.show("Results");

        IJ.log(String.format(
            Locale.US,
            "NPS ROI Batch: image=%s rois=%d slices=%d ROI size=%dx%d pixel=%.9f mm df=%.9f cy/mm Nyquist=%.9f cy/mm trend=%s source=%s",
            imp.getTitle(),
            rois.size(),
            result.perSlice.length,
            rois.get(0).roi.getBounds().width,
            rois.get(0).roi.getBounds().height,
            params.pixelSizeMm,
            result.frequencyStep,
            result.nyquist,
            params.trendRemoval,
            params.geometrySource
        ));
    }

    private void showMeasuredRoiPreview(ImagePlus imp, List<RoiEntry> rois) {
        Overlay overlay = new Overlay();
        for (int i = 0; i < rois.size(); i++) {
            RoiEntry entry = rois.get(i);
            try {
                Rectangle bounds = entry.roi.getBounds();
                Roi measured = new Roi(bounds.x, bounds.y, bounds.width, bounds.height);
                measured.setName("NPS measured ROI: " + entry.displayName);
                measured.setStrokeColor(roiColor(i, rois.size()));
                measured.setStrokeWidth(2.0);
                measured.setPosition(0);
                overlay.add(measured);
            } catch (RuntimeException e) {
                IJ.log("NPS ROI Batch: ROI overlay skipped for " + entry.displayName + ": " + e.getMessage());
            }
        }
        imp.setOverlay(overlay);
        imp.updateAndDraw();
        IJ.log("NPS ROI Batch ROI preview shown on the image. Colored rectangles mark all measured ROIs and are visible on every stack slice.");
    }

    private Color roiColor(int index, int total) {
        if (total <= 1) {
            return Color.YELLOW;
        }
        float hue = (float) ((index * 0.61803398875) % 1.0);
        return Color.getHSBColor(hue, 0.84f, 0.95f);
    }

    private void showNpsPlot(NpsResult result) {
        int start = firstOutputBin(result);
        int count = result.frequencies.length - start;
        if (count <= 0) {
            IJ.log("NPS ROI Batch plot skipped: no non-zero frequency bins available.");
            return;
        }

        double[] yRange = positiveNpsRange(result, start);
        if (yRange == null) {
            IJ.log("NPS ROI Batch plot skipped: no positive non-zero frequency NPS values available for logarithmic display.");
            return;
        }

        String title = "NPS ROI Batch";
        closeWindowByTitle(title);
        Plot plot = new Plot(title, "Spatial frequency (cycles/mm)", "NPS (HU^2 mm^2)");
        plot.setFrameSize(760, 440);
        plot.setAxisYLog(true);
        plot.setLimits(0.0, result.nyquist, yRange[0], yRange[1]);
        plot.setLineWidth(1);
        for (int s = 0; s < result.perSlice.length; s++) {
            plot.setColor(sliceColor(s, result.perSlice.length));
            drawSolidCurve(plot, result.frequencies, result.perSlice[s], start);
        }
        plot.setLineWidth(3);
        plot.setColor(Color.BLACK);
        drawSolidCurve(plot, result.frequencies, result.average, start);
        plot.setLineWidth(1);
        PlotWindow window = plot.show();
        if (window != null) {
            window.toFront();
        }
        IJ.log("NPS ROI Batch plot shown: thin colored curves are ROI-averaged per-slice NPS; thick black curve is the slice-average NPS. The 0 cycles/mm DC bin is omitted from data and Excel output, while the x-axis starts at 0 cycles/mm.");
    }

    private double[] positiveNpsRange(NpsResult result, int start) {
        double min = Double.POSITIVE_INFINITY;
        double max = Double.NEGATIVE_INFINITY;
        for (int s = 0; s < result.perSlice.length; s++) {
            for (int r = start; r < result.perSlice[s].length; r++) {
                double value = result.perSlice[s][r];
                if (isPositiveFinite(value)) {
                    if (value < min) min = value;
                    if (value > max) max = value;
                }
            }
        }
        for (int r = start; r < result.average.length; r++) {
            double value = result.average[r];
            if (isPositiveFinite(value)) {
                if (value < min) min = value;
                if (value > max) max = value;
            }
        }
        if (!Double.isFinite(min) || !Double.isFinite(max)) {
            return null;
        }
        if (min == max) {
            min /= 10.0;
            max *= 10.0;
        } else {
            min = Math.pow(10.0, Math.floor(Math.log10(min)));
            max = Math.pow(10.0, Math.ceil(Math.log10(max)));
        }
        if (min <= 0.0) {
            min = max / 1.0e6;
        }
        return new double[]{min, max};
    }

    private boolean isPositiveFinite(double value) {
        return value > 0.0 && !Double.isNaN(value) && !Double.isInfinite(value);
    }

    private Color sliceColor(int index, int total) {
        if (total <= 1) {
            return new Color(80, 120, 200);
        }
        float hue = (float) ((index * 0.61803398875) % 1.0);
        return Color.getHSBColor(hue, 0.72f, 0.78f);
    }

    private void drawSolidCurve(Plot plot, double[] x, double[] y, int start) {
        for (int r = start; r < x.length - 1 && r < y.length - 1; r++) {
            double x0 = x[r];
            double y0 = y[r];
            double x1 = x[r + 1];
            double y1 = y[r + 1];
            if (!isPositiveFinite(y0) || !isPositiveFinite(y1)) {
                continue;
            }
            plot.drawLine(x0, y0, x1, y1);
        }
    }

    private void closeWindowByTitle(String title) {
        Window existing = WindowManager.getWindow(title);
        if (existing != null) {
            existing.dispose();
        }
    }

    private boolean exportWorkbook(BatchParameters params, NpsResult result) {
        SaveDialog sd = new SaveDialog(
            "Save NPS ROI Batch Excel workbook",
            "NPS_2DFFT_ROIavg" + result.roiCount + "_slices" + params.startSlice + "-" + params.endSlice,
            ".xlsx"
        );
        String fileName = sd.getFileName();
        String directory = sd.getDirectory();
        if (fileName == null || directory == null) {
            IJ.log("NPS ROI Batch Excel export cancelled.");
            return false;
        }

        Path workbookPath = Paths.get(directory, ensureXlsxExtension(fileName));
        try {
            writeXlsxWorkbook(workbookPath, result);
            IJ.log("NPS ROI Batch Excel workbook saved to " + workbookPath.toString());
            return true;
        } catch (IOException e) {
            IJ.error("Export failed", e.getMessage());
            return false;
        }
    }

    private void writeXlsxWorkbook(Path path, NpsResult result) throws IOException {
        ZipOutputStream zip = new ZipOutputStream(Files.newOutputStream(path));
        try {
            putZipEntry(zip, "[Content_Types].xml", contentTypesXml());
            putZipEntry(zip, "_rels/.rels", rootRelationshipsXml());
            putZipEntry(zip, "xl/workbook.xml", workbookXml());
            putZipEntry(zip, "xl/_rels/workbook.xml.rels", workbookRelationshipsXml());
            putZipEntry(zip, "xl/worksheets/sheet1.xml", npsSheetXml(result));
        } finally {
            zip.close();
        }
    }

    private void putZipEntry(ZipOutputStream zip, String name, String text) throws IOException {
        zip.putNextEntry(new ZipEntry(name));
        byte[] bytes = text.getBytes(StandardCharsets.UTF_8);
        zip.write(bytes, 0, bytes.length);
        zip.closeEntry();
    }

    private String contentTypesXml() {
        return "<?xml version=\"1.0\" encoding=\"UTF-8\"?>"
            + "<Types xmlns=\"http://schemas.openxmlformats.org/package/2006/content-types\">"
            + "<Default Extension=\"rels\" ContentType=\"application/vnd.openxmlformats-package.relationships+xml\"/>"
            + "<Default Extension=\"xml\" ContentType=\"application/xml\"/>"
            + "<Override PartName=\"/xl/workbook.xml\" ContentType=\"application/vnd.openxmlformats-officedocument.spreadsheetml.sheet.main+xml\"/>"
            + "<Override PartName=\"/xl/worksheets/sheet1.xml\" ContentType=\"application/vnd.openxmlformats-officedocument.spreadsheetml.worksheet+xml\"/>"
            + "</Types>";
    }

    private String rootRelationshipsXml() {
        return "<?xml version=\"1.0\" encoding=\"UTF-8\"?>"
            + "<Relationships xmlns=\"http://schemas.openxmlformats.org/package/2006/relationships\">"
            + "<Relationship Id=\"rId1\" Type=\"http://schemas.openxmlformats.org/officeDocument/2006/relationships/officeDocument\" Target=\"xl/workbook.xml\"/>"
            + "</Relationships>";
    }

    private String workbookXml() {
        return "<?xml version=\"1.0\" encoding=\"UTF-8\"?>"
            + "<workbook xmlns=\"http://schemas.openxmlformats.org/spreadsheetml/2006/main\" "
            + "xmlns:r=\"http://schemas.openxmlformats.org/officeDocument/2006/relationships\">"
            + "<sheets>"
            + "<sheet name=\"NPS\" sheetId=\"1\" r:id=\"rId1\"/>"
            + "</sheets>"
            + "</workbook>";
    }

    private String workbookRelationshipsXml() {
        return "<?xml version=\"1.0\" encoding=\"UTF-8\"?>"
            + "<Relationships xmlns=\"http://schemas.openxmlformats.org/package/2006/relationships\">"
            + "<Relationship Id=\"rId1\" Type=\"http://schemas.openxmlformats.org/officeDocument/2006/relationships/worksheet\" Target=\"worksheets/sheet1.xml\"/>"
            + "</Relationships>";
    }

    private String npsSheetXml(NpsResult result) {
        int start = firstOutputBin(result);
        StringBuilder sb = new StringBuilder();
        sb.append("<?xml version=\"1.0\" encoding=\"UTF-8\"?>")
          .append("<worksheet xmlns=\"http://schemas.openxmlformats.org/spreadsheetml/2006/main\">")
          .append("<sheetData>");

        appendNpsHeaderRow(sb, result.sliceNumbers);
        int row = 2;
        for (int r = start; r < result.frequencies.length; r++) {
            sb.append("<row r=\"").append(row).append("\">");
            appendNumberCell(sb, row, 1, result.frequencies[r]);
            appendNumberCell(sb, row, 2, result.average[r]);
            for (int s = 0; s < result.sliceNumbers.length; s++) {
                appendNumberCell(sb, row, s + 3, result.perSlice[s][r]);
            }
            sb.append("</row>");
            row++;
        }

        sb.append("</sheetData></worksheet>");
        return sb.toString();
    }

    private void appendNpsHeaderRow(StringBuilder sb, int[] sliceNumbers) {
        sb.append("<row r=\"1\">");
        appendInlineStringCell(sb, 1, 1, "Spatial_frequency_cycles_per_mm");
        appendInlineStringCell(sb, 1, 2, "Average_NPS_HU2_mm2");
        for (int s = 0; s < sliceNumbers.length; s++) {
            appendInlineStringCell(sb, 1, s + 3, "Slice_" + sliceNumbers[s] + "_NPS_HU2_mm2");
        }
        sb.append("</row>");
    }

    private void appendInlineStringCell(StringBuilder sb, int row, int col, String value) {
        sb.append("<c r=\"").append(cellRef(row, col)).append("\" t=\"inlineStr\"><is><t>")
          .append(xmlEscape(value))
          .append("</t></is></c>");
    }

    private void appendNumberCell(StringBuilder sb, int row, int col, double value) {
        sb.append("<c r=\"").append(cellRef(row, col)).append("\">");
        if (Double.isFinite(value)) {
            sb.append("<v>").append(numberCell(value)).append("</v>");
        }
        sb.append("</c>");
    }

    private String cellRef(int row, int col) {
        return columnName(col) + row;
    }

    private String columnName(int col) {
        StringBuilder sb = new StringBuilder();
        int value = col;
        while (value > 0) {
            int rem = (value - 1) % 26;
            sb.insert(0, (char) ('A' + rem));
            value = (value - 1) / 26;
        }
        return sb.toString();
    }

    private String xmlEscape(String value) {
        if (value == null) {
            return "";
        }
        StringBuilder sb = new StringBuilder(value.length());
        for (int i = 0; i < value.length(); i++) {
            char ch = value.charAt(i);
            if (ch == '&') {
                sb.append("&amp;");
            } else if (ch == '<') {
                sb.append("&lt;");
            } else if (ch == '>') {
                sb.append("&gt;");
            } else if (ch == '"') {
                sb.append("&quot;");
            } else if (ch >= 0x20 || ch == '\t' || ch == '\n' || ch == '\r') {
                sb.append(ch);
            }
        }
        return sb.toString();
    }

    private String numberCell(double value) {
        return Double.isNaN(value) || Double.isInfinite(value)
            ? ""
            : String.format(Locale.US, "%.12g", value);
    }

    private String ensureXlsxExtension(String fileName) {
        return fileName.toLowerCase(Locale.US).endsWith(".xlsx") ? fileName : fileName + ".xlsx";
    }

    private int firstOutputBin(NpsResult result) {
        return result != null && result.frequencies.length > FIRST_OUTPUT_BIN ? FIRST_OUTPUT_BIN : 0;
    }

    private int clampSlice(int slice, int stackSize) {
        return Math.max(1, Math.min(stackSize, slice));
    }

    private int[] sliceNumbers(int startSlice, int endSlice, int step) {
        int count = ((endSlice - startSlice) / step) + 1;
        int[] slices = new int[count];
        int value = startSlice;
        for (int i = 0; i < count; i++) {
            slices[i] = value;
            value += step;
        }
        return slices;
    }

    private String baseName(String entryName) {
        String normalized = entryName.replace('\\', '/');
        int slash = normalized.lastIndexOf('/');
        String name = slash >= 0 ? normalized.substring(slash + 1) : normalized;
        if (name.toLowerCase(Locale.US).endsWith(".roi")) {
            name = name.substring(0, name.length() - 4);
        }
        return name;
    }

    private void clearStatusAndProgress() {
        IJ.showProgress(1.0);
        IJ.showStatus("");
    }

    private Geometry inferGeometry(ImagePlus imp, Calibration cal, int width, int height) {
        String info = dicomInfo(imp);
        // Pixel Spacing describes the sampled image grid. Header Rows/Columns may
        // be stale after a crop, so use the current image width for its FOV.
        String spacingPayload = findDicomPayload(info, "0028,0030", "Pixel Spacing");
        if (spacingPayload != null) {
            double[] spacing = parseDicomDoubles(info, "0028,0030", "Pixel Spacing");
            if (spacing.length != 2 || !isPositive(Double.valueOf(spacing[0]))
                || !isPositive(Double.valueOf(spacing[1]))) {
                throw new IllegalArgumentException("DICOM Pixel Spacing must contain two finite positive values.");
            }
            if (!nearlyEqual(spacing[0], spacing[1])) {
                throw new IllegalArgumentException("Non-square pixels are not supported by radial NPS. DICOM Pixel Spacing differs between rows and columns.");
            }
            double pixel = spacing[1];
            return new Geometry(pixel, pixel * width, width, true,
                String.format(Locale.US, "DICOM Pixel Spacing (0028,0030)=%.9g mm; current image width=%d", pixel, width));
        }

        Double reconstructionDiameter = parseDicomDouble(info, "0018,1100", "Reconstruction Diameter");
        if (isPositive(reconstructionDiameter)) {
            Integer columns = parseDicomInteger(info, "0028,0011", "Columns");
            Integer rows = parseDicomInteger(info, "0028,0010", "Rows");
            if ((isPositive(columns) && columns.intValue() != width)
                || (isPositive(rows) && rows.intValue() != height)) {
                throw new IllegalArgumentException("DICOM dimensions differ from the current image. Pixel Spacing is required to determine its pixel size.");
            }
            if (width != height) {
                throw new IllegalArgumentException("Reconstruction Diameter alone cannot determine square-pixel geometry for a non-square image.");
            }
            double fovMm = reconstructionDiameter.doubleValue();
            return new Geometry(fovMm / width, fovMm, width, true,
                String.format(Locale.US, "DICOM Reconstruction Diameter (0018,1100)=%.9g mm / current image width=%d", fovMm, width));
        }

        double pixel = inferCalibrationPixelSizeMm(cal);
        return new Geometry(pixel, pixel * width, width, false,
            String.format(Locale.US, "DICOM geometry not found; calibration/default estimate %.9g mm/pixel", pixel));
    }

    private double inferCalibrationPixelSizeMm(Calibration cal) {
        if (cal == null) {
            return 1.0;
        }
        if (!isPositive(Double.valueOf(cal.pixelWidth)) || !isPositive(Double.valueOf(cal.pixelHeight))) {
            throw new IllegalArgumentException("Image calibration must have finite positive pixel width and height.");
        }
        if (!nearlyEqual(cal.pixelWidth, cal.pixelHeight)) {
            throw new IllegalArgumentException("Non-square calibrated pixels are not supported by radial NPS.");
        }
        String unit = cal.getUnit();
        if (unit == null || unit.equalsIgnoreCase("pixel") || unit.equalsIgnoreCase("pixels")) {
            return 1.0; // Uncalibrated input: the dialog explicitly requests geometry confirmation.
        }
        double factor;
        if (unit.equalsIgnoreCase("mm") || unit.equalsIgnoreCase("millimeter") || unit.equalsIgnoreCase("millimeters")) {
            factor = 1.0;
        } else if (unit.equalsIgnoreCase("cm")) {
            factor = 10.0;
        } else if (unit.equalsIgnoreCase("um") || unit.equalsIgnoreCase("micron") || unit.equalsIgnoreCase("microns")
            || unit.equalsIgnoreCase("micrometer") || unit.equalsIgnoreCase("micrometers") || unit.equals("\u00b5m") || unit.equals("\u03bcm")) {
            factor = 0.001;
        } else {
            throw new IllegalArgumentException("Unsupported spatial calibration unit: " + unit + ". Calibrate the image in mm before measuring NPS.");
        }
        double pixel = cal.pixelWidth * factor;
        if (!isPositive(Double.valueOf(pixel))) {
            throw new IllegalArgumentException("Calibrated pixel size in mm is invalid.");
        }
        return pixel;
    }

    private String dicomInfo(ImagePlus imp) {
        if (imp == null) {
            return "";
        }
        StringBuilder sb = new StringBuilder();
        appendInfo(sb, imp.getInfoProperty());
        Object property = imp.getProperty("Info");
        if (property != null) {
            appendInfo(sb, property.toString());
        }
        FileInfo fileInfo = imp.getOriginalFileInfo();
        if (fileInfo != null) {
            appendInfo(sb, fileInfo.info);
        }
        if (imp.getStack() != null && imp.getStackSize() > 0) {
            int current = Math.max(1, Math.min(imp.getCurrentSlice(), imp.getStackSize()));
            appendInfo(sb, imp.getStack().getSliceLabel(current));
            if (current != 1) {
                appendInfo(sb, imp.getStack().getSliceLabel(1));
            }
        }
        return sb.toString();
    }

    private void appendInfo(StringBuilder sb, String text) {
        if (text == null || text.trim().length() == 0) {
            return;
        }
        if (sb.length() > 0) {
            sb.append('\n');
        }
        sb.append(text);
    }

    private Integer parseDicomInteger(String info, String tag, String label) {
        Double value = parseDicomDouble(info, tag, label);
        if (!isPositive(value)) {
            return null;
        }
        return Integer.valueOf((int) Math.round(value.doubleValue()));
    }

    private Double parseDicomDouble(String info, String tag, String label) {
        double[] values = parseDicomDoubles(info, tag, label);
        return values.length == 0 ? null : Double.valueOf(values[0]);
    }

    private double[] parseDicomDoubles(String info, String tag, String label) {
        String payload = findDicomPayload(info, tag, label);
        if (payload == null) {
            return new double[0];
        }
        List<Double> values = new ArrayList<>();
        Matcher matcher = NUMBER_PATTERN.matcher(payload);
        while (matcher.find()) {
            try {
                values.add(Double.valueOf(Double.parseDouble(matcher.group())));
            } catch (NumberFormatException ignored) {
                // Ignore malformed numeric fragments in vendor-specific DICOM text.
            }
        }
        double[] parsed = new double[values.size()];
        for (int i = 0; i < values.size(); i++) {
            parsed[i] = values.get(i).doubleValue();
        }
        return parsed;
    }

    private String findDicomPayload(String info, String tag, String label) {
        if (info == null || info.length() == 0) {
            return null;
        }
        String normalizedTag = tag.toLowerCase(Locale.US);
        String tagWithoutComma = normalizedTag.replace(",", "");
        String[] lines = info.split("\\r?\\n");
        for (String line : lines) {
            String normalizedLine = line.toLowerCase(Locale.US).replace("(", "").replace(")", "").replace(" ", "");
            int tagIndex = normalizedLine.indexOf(normalizedTag);
            int tagLength = normalizedTag.length();
            if (tagIndex < 0) {
                tagIndex = normalizedLine.indexOf(tagWithoutComma);
                tagLength = tagWithoutComma.length();
            }
            if (tagIndex >= 0) {
                int rawStart = approximateRawTagEnd(line, tag, tagIndex, tagLength);
                return line.substring(Math.min(rawStart, line.length()));
            }
        }

        String normalizedLabel = label == null ? "" : label.toLowerCase(Locale.US);
        if (normalizedLabel.length() == 0) {
            return null;
        }
        if (normalizedLabel.equals("pixel spacing")) {
            // Do not confuse detector/imager spacing or a similarly named tag
            // with the reconstructed image's Pixel Spacing (0028,0030).
            Pattern exactLabel = Pattern.compile("^\\s*Pixel\\s+Spacing\\s*[:=]\\s*(.*)$", Pattern.CASE_INSENSITIVE);
            for (String line : lines) {
                Matcher match = exactLabel.matcher(line);
                if (match.matches()) {
                    return match.group(1);
                }
            }
            return null;
        }
        for (String line : lines) {
            String lower = line.toLowerCase(Locale.US);
            int labelIndex = lower.indexOf(normalizedLabel);
            if (labelIndex >= 0) {
                return line.substring(labelIndex + normalizedLabel.length());
            }
        }
        return null;
    }

    private int approximateRawTagEnd(String line, String tag, int normalizedTagIndex, int normalizedTagLength) {
        int rawStart = 0;
        int compactCount = 0;
        int target = normalizedTagIndex + normalizedTagLength;
        for (int i = 0; i < line.length(); i++) {
            char c = line.charAt(i);
            if (c == '(' || c == ')' || c == ' ') {
                continue;
            }
            compactCount++;
            if (compactCount >= target) {
                rawStart = i + 1;
                break;
            }
        }
        if (rawStart <= 0) {
            int directIndex = line.toLowerCase(Locale.US).indexOf(tag.toLowerCase(Locale.US));
            rawStart = directIndex >= 0 ? directIndex + tag.length() : line.length();
        }
        while (rawStart < line.length()) {
            char c = line.charAt(rawStart);
            if (c == ')' || c == ':' || c == '=' || Character.isWhitespace(c)) {
                rawStart++;
            } else {
                break;
            }
        }
        return rawStart;
    }

    private boolean isPowerOfTwo(int value) {
        return value > 0 && (value & (value - 1)) == 0;
    }

    private boolean isPositive(Number value) {
        return value != null && value.doubleValue() > 0.0 && !Double.isNaN(value.doubleValue()) && !Double.isInfinite(value.doubleValue());
    }

    private boolean nearlyEqual(double a, double b) {
        double tolerance = Math.max(1e-6, Math.max(Math.abs(a), Math.abs(b)) * 1e-4);
        return Math.abs(a - b) <= tolerance;
    }

    private static class Geometry {
        final double pixelSizeMm;
        final double fovMm;
        final int matrix;
        final boolean isDicomDerived;
        final String source;

        Geometry(double pixelSizeMm, double fovMm, int matrix, boolean isDicomDerived, String source) {
            this.pixelSizeMm = pixelSizeMm;
            this.fovMm = fovMm;
            this.matrix = matrix;
            this.isDicomDerived = isDicomDerived;
            this.source = source;
        }
    }

    private static class RoiEntry {
        final Roi roi;
        final String name;
        final String entryName;
        String displayName;

        RoiEntry(Roi roi, String name, String entryName) {
            this.roi = roi;
            this.name = name;
            this.entryName = entryName;
            this.displayName = name;
        }
    }

    private static class BatchParameters {
        final double pixelSizeMm;
        final double fovMm;
        final int matrix;
        final boolean geometryDicomDerived;
        final String geometrySource;
        final int[] sliceNumbers;
        final int startSlice;
        final int endSlice;
        final int sliceStep;
        final String trendRemoval;
        final boolean useHannWindow;

        BatchParameters(double pixelSizeMm, double fovMm, int matrix, boolean geometryDicomDerived,
                        String geometrySource, int[] sliceNumbers, int startSlice, int endSlice,
                        int sliceStep, String trendRemoval, boolean useHannWindow) {
            this.pixelSizeMm = pixelSizeMm;
            this.fovMm = fovMm;
            this.matrix = matrix;
            this.geometryDicomDerived = geometryDicomDerived;
            this.geometrySource = geometrySource;
            this.sliceNumbers = sliceNumbers;
            this.startSlice = startSlice;
            this.endSlice = endSlice;
            this.sliceStep = sliceStep;
            this.trendRemoval = trendRemoval;
            this.useHannWindow = useHannWindow;
        }

        NpsParameters forRoi(RoiEntry roi) {
            return new NpsParameters(
                pixelSizeMm,
                fovMm,
                matrix,
                geometryDicomDerived,
                geometrySource,
                roi.roi.getBounds(),
                sliceNumbers,
                trendRemoval,
                useHannWindow
            );
        }
    }

    private static class NpsParameters {
        final double pixelSizeMm;
        final double fovMm;
        final int matrix;
        final boolean geometryDicomDerived;
        final String geometrySource;
        final Rectangle roiBounds;
        final int[] sliceNumbers;
        final String trendRemoval;
        final boolean useHannWindow;

        NpsParameters(double pixelSizeMm, double fovMm, int matrix, boolean geometryDicomDerived,
                      String geometrySource, Rectangle roiBounds, int[] sliceNumbers,
                      String trendRemoval, boolean useHannWindow) {
            this.pixelSizeMm = pixelSizeMm;
            this.fovMm = fovMm;
            this.matrix = matrix;
            this.geometryDicomDerived = geometryDicomDerived;
            this.geometrySource = geometrySource;
            this.roiBounds = roiBounds;
            this.sliceNumbers = sliceNumbers;
            this.trendRemoval = trendRemoval;
            this.useHannWindow = useHannWindow;
        }
    }

    private static class NpsResult {
        final double[] frequencies;
        final double[][] perSlice;
        final double[] average;
        final int[] sliceNumbers;
        final double frequencyStep;
        final double nyquist;
        final NpsParameters parameters;
        final int roiCount;

        NpsResult(double[] frequencies, double[][] perSlice, double[] average, int[] sliceNumbers,
                  double frequencyStep, double nyquist, NpsParameters parameters, int roiCount) {
            this.frequencies = frequencies;
            this.perSlice = perSlice;
            this.average = average;
            this.sliceNumbers = sliceNumbers;
            this.frequencyStep = frequencyStep;
            this.nyquist = nyquist;
            this.parameters = parameters;
            this.roiCount = roiCount;
        }
    }
}
