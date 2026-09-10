import ij.IJ;
import ij.ImagePlus;
import ij.WindowManager;
import ij.gui.GenericDialog;
import ij.gui.ImageCanvas;
import ij.gui.ImageRoi;
import ij.gui.NonBlockingGenericDialog;
import ij.gui.Overlay;
import ij.gui.Plot;
import ij.gui.PlotWindow;
import ij.gui.Roi;
import ij.gui.YesNoCancelDialog;
import ij.io.FileInfo;
import ij.io.SaveDialog;
import ij.measure.Calibration;
import ij.measure.ResultsTable;
import ij.plugin.PlugIn;
import ij.process.FHT;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;

import java.awt.Choice;
import java.awt.Color;
import java.awt.Rectangle;
import java.awt.TextField;
import java.awt.Window;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.awt.event.TextEvent;
import java.awt.event.TextListener;
import java.util.ArrayList;
import java.util.List;
import java.util.Locale;
import java.util.Vector;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

/**
 * Noise power spectrum measurement with a calibrated 2D FFT method.
 */
public class NPS_2D_FFT implements PlugIn {
    private static final String VERSION = "2026-06-19-nps-quadratic-detrend";
    private static final int DEFAULT_ROI_SIZE = 128;
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
        IJ.log("NPS Calculator " + VERSION + " loaded as class " + getClass().getName());
        ImagePlus imp = IJ.getImage();
        if (imp == null) {
            IJ.error("No image", "Open an image first.");
            return;
        }
        if (imp.getType() == ImagePlus.COLOR_RGB) {
            IJ.error("Unsupported image", "NPS should be calculated from a grayscale CT image.");
            return;
        }

        Calibration cal = imp.getCalibration();
        Geometry geometry = inferGeometry(imp, cal, imp.getWidth(), imp.getHeight());
        while (true) {
            NpsParameters params = readParameters(imp, geometry);
            if (params == null) {
                return;
            }

            NpsResult result = computeNps(imp, cal, params);
            showSummary(imp, params, result);
            showNpsPlot(result);
            boolean csvSaved = exportCsv(result);
            if (confirmFinishAfterCalculation("NPS", csvSaved)) {
                return;
            }
        }
    }

    private NpsParameters readParameters(ImagePlus imp, Geometry geometry) {
        int width = imp.getWidth();
        int height = imp.getHeight();
        int matrixLimit = Math.max(1, geometry.matrix);
        int defaultRoi = largestPowerOfTwoAtMost(Math.min(Math.min(Math.min(width, height), matrixLimit), DEFAULT_ROI_SIZE));
        if (defaultRoi < 8) {
            IJ.error("Image too small", "Image must support at least an 8x8 square ROI.");
            return null;
        }

        boolean hasUsableRoi = isUsableSquarePowerOfTwoRoi(imp.getRoi(), width, height, matrixLimit);
        if (hasUsableRoi) {
            defaultRoi = imp.getRoi().getBounds().width;
        }
        Rectangle defaultBounds = defaultRoiBounds(imp, defaultRoi, hasUsableRoi);
        double defaultCenterX = defaultBounds.x + defaultBounds.width / 2.0;
        double defaultCenterY = defaultBounds.y + defaultBounds.height / 2.0;
        String[] roiSizeChoices = allowedPowerOfTwoRoiSizes(width, height, matrixLimit);
        GenericDialog gd = new NonBlockingGenericDialog("NPS 2D FFT Parameters");
        gd.addMessage(String.format(
            Locale.US,
            "Geometry: FOV %.9g mm / matrix %d = %.9g mm/pixel%nSource: %s",
            geometry.fovMm,
            geometry.matrix,
            geometry.pixelSizeMm,
            geometry.source
        ));
        gd.addMessage(String.format(
            Locale.US,
            "Frequency step uses DICOM-derived pixel size and the FFT ROI length: df = 1 / (ROI size * pixel size).%nFor %d px ROI: df = %.9g cycles/mm.",
            defaultRoi,
            1.0 / (defaultRoi * geometry.pixelSizeMm)
        ));
        if (!geometry.isDicomDerived) {
            gd.addMessage("DICOM FOV/pixel-spacing tags were not found. Confirm or edit FOV and matrix below.");
        }
        gd.addNumericField("FOV (mm)", geometry.fovMm, 6);
        gd.addNumericField("Matrix size (pixels)", geometry.matrix, 0);
        gd.addChoice("Square ROI size (pixels; power of 2, <= matrix)", roiSizeChoices, Integer.toString(defaultRoi));
        gd.addNumericField("ROI center X (pixels)", defaultCenterX, 2);
        gd.addNumericField("ROI center Y (pixels)", defaultCenterY, 2);
        gd.addMessage("Move the cursor over the image to read x/y/value. Click the image to set ROI center X/Y.");
        gd.addMessage("Current ROI option accepts only a rectangular square ROI whose side length is a power of two; oval/circular ROIs are rejected.");
        gd.addCheckbox("Use current rectangular square ROI instead of center X/Y", hasUsableRoi);
        gd.addChoice("Trend removal before FFT", TREND_REMOVAL_CHOICES, TREND_REMOVAL_QUADRATIC);
        gd.addCheckbox("Apply 2D Hann window", false);
        RoiPreviewController preview = installRoiPreview(imp, gd);
        ImageCoordinateProbe probe = installImageCoordinateProbe(imp, gd, imp.getCalibration(), preview);
        gd.showDialog();
        if (probe != null) {
            probe.uninstall();
        }
        if (preview != null) {
            preview.uninstall();
        }
        if (gd.wasCanceled()) {
            return null;
        }

        double enteredFov = gd.getNextNumber();
        double enteredMatrix = gd.getNextNumber();
        double enteredCenterX = gd.getNextNumber();
        double enteredCenterY = gd.getNextNumber();
        String roiSizeChoice = gd.getNextChoice();
        boolean useCurrentRoi = gd.getNextBoolean();
        String trendRemoval = gd.getNextChoice();
        boolean useHann = gd.getNextBoolean();
        if (gd.invalidNumber() || enteredFov <= 0 || enteredMatrix <= 0
            || Double.isNaN(enteredFov) || Double.isNaN(enteredMatrix)
            || Double.isNaN(enteredCenterX) || Double.isNaN(enteredCenterY)) {
            IJ.error("Invalid input", "FOV, matrix, and ROI center must be valid numbers.");
            return null;
        }

        int matrix = Math.max(1, (int) Math.round(enteredMatrix));
        int roiSize = Integer.parseInt(roiSizeChoice);
        Rectangle roiBounds;
        if (useCurrentRoi) {
            Roi roi = imp.getRoi();
            if (!isUsableSquarePowerOfTwoRoi(roi, width, height, matrix)) {
                IJ.error("Invalid ROI", "Current ROI must be a rectangular square, power-of-two sized, and inside the image. Oval/circular ROIs are not valid for this 2D FFT method.");
                return null;
            }
            roiBounds = roi.getBounds();
            roiSize = roiBounds.width;
        } else {
            if (roiSize > matrix) {
                IJ.error("Invalid ROI size", "ROI size must not exceed the matrix size.");
                return null;
            }
            if (roiSize > width || roiSize > height) {
                IJ.error("Invalid ROI size", "ROI size must fit inside the image.");
                return null;
            }
            int x0 = (int) Math.round(enteredCenterX - roiSize / 2.0);
            int y0 = (int) Math.round(enteredCenterY - roiSize / 2.0);
            roiBounds = new Rectangle(x0, y0, roiSize, roiSize);
            if (roiBounds.x < 0 || roiBounds.y < 0 || roiBounds.x + roiSize > width || roiBounds.y + roiSize > height) {
                IJ.error("Invalid ROI position", "The selected ROI center and size place the ROI outside the image.");
                return null;
            }
        }

        double pixelSizeMm = enteredFov / matrix;
        boolean unchangedGeometry = nearlyEqual(enteredFov, geometry.fovMm) && matrix == geometry.matrix;
        String source = unchangedGeometry
            ? geometry.source
            : String.format(Locale.US, "Manual dialog input: FOV=%.9g mm, matrix=%d", enteredFov, matrix);
        return new NpsParameters(
            pixelSizeMm,
            enteredFov,
            matrix,
            geometry.isDicomDerived && unchangedGeometry,
            source,
            roiBounds,
            trendRemoval,
            useHann
        );
    }

    private Rectangle defaultRoiBounds(ImagePlus imp, int defaultRoi, boolean hasUsableRoi) {
        if (hasUsableRoi) {
            return imp.getRoi().getBounds();
        }
        int x0 = (imp.getWidth() - defaultRoi) / 2;
        int y0 = (imp.getHeight() - defaultRoi) / 2;
        return new Rectangle(x0, y0, defaultRoi, defaultRoi);
    }

    private RoiPreviewController installRoiPreview(ImagePlus imp, GenericDialog gd) {
        Vector<?> fields = gd.getNumericFields();
        Vector<?> choices = gd.getChoices();
        if (fields == null || fields.size() < 4 || choices == null || choices.size() < 1) {
            return null;
        }
        Choice roiSizeChoice = (Choice) choices.get(0);
        TextField matrixField = (TextField) fields.get(1);
        TextField centerXField = (TextField) fields.get(2);
        TextField centerYField = (TextField) fields.get(3);
        RoiPreviewController preview = new RoiPreviewController(imp, roiSizeChoice, matrixField, centerXField, centerYField);
        preview.install();
        return preview;
    }

    private ImageCoordinateProbe installImageCoordinateProbe(ImagePlus imp, GenericDialog gd, Calibration cal, RoiPreviewController preview) {
        ImageCanvas canvas = imp.getCanvas();
        if (canvas == null) {
            return null;
        }
        Vector<?> fields = gd.getNumericFields();
        if (fields == null || fields.size() < 4) {
            return null;
        }
        TextField centerXField = (TextField) fields.get(2);
        TextField centerYField = (TextField) fields.get(3);
        ImageCoordinateProbe probe = new ImageCoordinateProbe(imp, canvas, cal, centerXField, centerYField, preview);
        canvas.addMouseListener(probe);
        canvas.addMouseMotionListener(probe);
        return probe;
    }

    private NpsResult computeNps(ImagePlus imp, Calibration cal, NpsParameters params) {
        int nSlices = imp.getStackSize();
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
        double[][] perSlice = new double[nSlices][maxRadius + 1];
        double[] average = new double[maxRadius + 1];

        for (int slice = 1; slice <= nSlices; slice++) {
            FloatProcessor fp = calibratedRoiToFloat(imp.getStack().getProcessor(slice), cal, params.roiBounds);
            applyTrendRemoval(fp, params.trendRemoval);
            applySeparableWindow(fp, window);

            FHT fht = new FHT(fp);
            fht.transform();
            FloatProcessor power = fht.getRawPowerSpectrum();
            fht.swapQuadrants(power);

            double[] radial = radialAverage(power, maxRadius);
            for (int r = 0; r <= maxRadius; r++) {
                double nps = radial[r] * pixelAreaMm2 / (roiSize * roiSize) / meanWindowSquared;
                perSlice[slice - 1][r] = nps;
                average[r] += nps;
            }
        }

        for (int r = 0; r <= maxRadius; r++) {
            average[r] /= nSlices;
        }
        return new NpsResult(frequencies, perSlice, average, frequencyStep, 1.0 / (2.0 * params.pixelSizeMm), params);
    }

    private FloatProcessor calibratedRoiToFloat(ImageProcessor ip, Calibration cal, Rectangle roi) {
        FloatProcessor fp = new FloatProcessor(roi.width, roi.height);
        for (int y = 0; y < roi.height; y++) {
            for (int x = 0; x < roi.width; x++) {
                double raw = ip.getf(roi.x + x, roi.y + y);
                double value = cal != null ? cal.getCValue(raw) : raw;
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
            IJ.log("NPS: 2D second-order trend fit failed; falling back to ROI mean subtraction.");
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

    private void showSummary(ImagePlus imp, NpsParameters params, NpsResult result) {
        ResultsTable rt = ResultsTable.getResultsTable();
        if (rt == null) {
            rt = new ResultsTable();
        }
        rt.incrementCounter();
        rt.addValue("Slices", result.perSlice.length);
        rt.addValue("PixelSize_mm", params.pixelSizeMm);
        rt.addValue("FOV_mm", params.fovMm);
        rt.addValue("Matrix", params.matrix);
        rt.addValue("Geometry_DICOM_Derived", params.geometryDicomDerived ? 1 : 0);
        rt.addValue("ROI_X", params.roiBounds.x);
        rt.addValue("ROI_Y", params.roiBounds.y);
        rt.addValue("ROI_Size_px", params.roiBounds.width);
        rt.addValue("Quadratic_Trend_Removed", TREND_REMOVAL_QUADRATIC.equals(params.trendRemoval) ? 1 : 0);
        rt.addValue("Mean_Subtracted", TREND_REMOVAL_NONE.equals(params.trendRemoval) ? 0 : 1);
        rt.addValue("Hann_Window_Applied", params.useHannWindow ? 1 : 0);
        rt.addValue("FFT_Length_mm", params.roiBounds.width * params.pixelSizeMm);
        rt.addValue("FrequencyStep_cycles_per_mm", result.frequencyStep);
        rt.addValue("Nyquist_cycles_per_mm", result.nyquist);
        rt.show("Results");

        IJ.log(String.format(
            Locale.US,
            "NPS: image=%s slices=%d ROI=%dx%d at (%d,%d) pixel=%.9f mm FOV=%.6f mm matrix=%d FFT length=%.9f mm df=1/(ROI*pixel)=%.9f cy/mm Nyquist=%.9f cy/mm trend=%s source=%s",
            imp.getTitle(),
            result.perSlice.length,
            params.roiBounds.width,
            params.roiBounds.height,
            params.roiBounds.x,
            params.roiBounds.y,
            params.pixelSizeMm,
            params.fovMm,
            params.matrix,
            params.roiBounds.width * params.pixelSizeMm,
            result.frequencyStep,
            result.nyquist,
            params.trendRemoval,
            params.geometrySource
        ));
    }

    private void showNpsPlot(NpsResult result) {
        int start = firstOutputBin(result);
        int count = result.frequencies.length - start;
        if (count <= 0) {
            IJ.log("NPS plot skipped: no non-zero frequency bins available.");
            return;
        }

        double[] yRange = positiveNpsRange(result, start);
        if (yRange == null) {
            IJ.log("NPS plot skipped: no positive non-zero frequency NPS values available for logarithmic display.");
            return;
        }

        String title = "NPS 2D FFT";
        closeWindowByTitle(title);
        Plot plot = new Plot(title, "Spatial frequency (cycles/mm)", "NPS (HU^2 mm^2)");
        plot.setFrameSize(720, 420);
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
        IJ.log("NPS plot shown with per-slice solid curves and a thick average curve. The 0 cycles/mm DC bin is omitted from data and CSV output, while the x-axis starts at 0 cycles/mm.");
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

    private boolean exportCsv(NpsResult result) {
        SaveDialog sd = new SaveDialog("Save NPS CSV", "NPS_2DFFT_results", ".csv");
        String dir = sd.getDirectory();
        String name = sd.getFileName();
        if (dir == null || name == null) {
            IJ.log("NPS CSV export cancelled.");
            clearStatusAndProgress();
            return false;
        }

        StringBuilder sb = new StringBuilder();
        sb.append("Frequency_cycles_per_mm");
        for (int s = 1; s <= result.perSlice.length; s++) {
            sb.append(",Slice").append(s).append("_NPS_HU2_mm2");
        }
        sb.append(",Average_NPS_HU2_mm2\n");
        for (int r = firstOutputBin(result); r < result.frequencies.length; r++) {
            sb.append(String.format(Locale.US, "%.9f", result.frequencies[r]));
            for (int s = 0; s < result.perSlice.length; s++) {
                sb.append(',').append(String.format(Locale.US, "%.12g", result.perSlice[s][r]));
            }
            sb.append(',').append(String.format(Locale.US, "%.12g", result.average[r])).append('\n');
        }
        IJ.saveString(sb.toString(), dir + name);
        IJ.log("NPS CSV saved: " + dir + name);
        clearStatusAndProgress();
        return true;
    }

    private int firstOutputBin(NpsResult result) {
        return result != null && result.frequencies.length > FIRST_OUTPUT_BIN ? FIRST_OUTPUT_BIN : 0;
    }

    private void closeWindowByTitle(String title) {
        Window existing = WindowManager.getWindow(title);
        if (existing != null) {
            existing.dispose();
        }
    }

    private boolean confirmFinishAfterCalculation(String analysisName, boolean csvSaved) {
        String csvStatus = csvSaved ? "CSVは保存済みです。" : "CSVは保存されていません。";
        YesNoCancelDialog dialog = new YesNoCancelDialog(
            IJ.getInstance(),
            analysisName + " calculation finished",
            "計算が完了しました。\n" + csvStatus + "\n" + analysisName + "をここで終了してよいですか？\n\n「続ける」またはキャンセルでは、次のROI/条件の入力に戻ります。",
            "終了",
            "続ける"
        );
        return dialog.yesPressed();
    }

    private void clearStatusAndProgress() {
        IJ.showProgress(1.0);
        IJ.showStatus("");
    }

    private Geometry inferGeometry(ImagePlus imp, Calibration cal, int width, int height) {
        String info = dicomInfo(imp);
        Integer rows = parseDicomInteger(info, "0028,0010", "Rows");
        Integer columns = parseDicomInteger(info, "0028,0011", "Columns");
        String matrixSource = "DICOM Rows/Columns (0028,0010/0028,0011)";
        if (!isPositive(rows)) {
            rows = Integer.valueOf(height);
            matrixSource = "ImageJ image dimensions";
        }
        if (!isPositive(columns)) {
            columns = Integer.valueOf(width);
            matrixSource = "ImageJ image dimensions";
        }

        Double reconstructionDiameter = parseDicomDouble(info, "0018,1100", "Reconstruction Diameter");
        if (isPositive(reconstructionDiameter)) {
            double fovMm = reconstructionDiameter.doubleValue();
            return new Geometry(
                fovMm / columns.doubleValue(),
                fovMm,
                columns.intValue(),
                true,
                String.format(Locale.US, "DICOM Reconstruction Diameter (0018,1100)=%.9g mm / matrix=%d (%s)", fovMm, columns, matrixSource)
            );
        }

        double[] pixelSpacing = parseDicomDoubles(info, "0028,0030", "Pixel Spacing");
        if (pixelSpacing.length > 0 && isPositive(Double.valueOf(pixelSpacing[0]))) {
            double rowSpacing = pixelSpacing[0];
            double columnSpacing = pixelSpacing.length > 1 && isPositive(Double.valueOf(pixelSpacing[1])) ? pixelSpacing[1] : rowSpacing;
            if (nearlyEqual(rowSpacing, columnSpacing)) {
                double fovMm = columnSpacing * columns.doubleValue();
                return new Geometry(
                    columnSpacing,
                    fovMm,
                    columns.intValue(),
                    true,
                    String.format(Locale.US, "DICOM Pixel Spacing (0028,0030)=%.9g mm * matrix=%d (%s)", columnSpacing, columns, matrixSource)
                );
            }
        }

        double calibrationPixel = inferCalibrationPixelSizeMm(cal);
        return new Geometry(
            calibrationPixel,
            calibrationPixel * width,
            width,
            false,
            String.format(Locale.US, "DICOM geometry not found; calibration/default estimate %.9g mm/pixel", calibrationPixel)
        );
    }

    private double inferCalibrationPixelSizeMm(Calibration cal) {
        if (cal == null) {
            return 1.0;
        }
        double candidate = cal.pixelWidth > 0 ? cal.pixelWidth : cal.pixelHeight;
        if (candidate <= 0) {
            return 1.0;
        }
        String unit = cal.getUnit();
        if (unit == null || unit.equalsIgnoreCase("pixel") || unit.equalsIgnoreCase("pixels")) {
            return 1.0;
        }
        if (unit.equalsIgnoreCase("cm")) {
            return candidate * 10.0;
        }
        if (unit.equalsIgnoreCase("um") || unit.equalsIgnoreCase("micron") || unit.equalsIgnoreCase("microns")) {
            return candidate / 1000.0;
        }
        return candidate;
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

    private boolean isUsableSquarePowerOfTwoRoi(Roi roi, int width, int height, int matrixLimit) {
        if (roi == null) {
            return false;
        }
        if (roi.getType() != Roi.RECTANGLE) {
            return false;
        }
        Rectangle bounds = roi.getBounds();
        return bounds.width == bounds.height
            && bounds.width > 0
            && bounds.width <= matrixLimit
            && isPowerOfTwo(bounds.width)
            && bounds.x >= 0
            && bounds.y >= 0
            && bounds.x + bounds.width <= width
            && bounds.y + bounds.height <= height;
    }

    private String[] allowedPowerOfTwoRoiSizes(int width, int height, int matrixLimit) {
        int max = largestPowerOfTwoAtMost(Math.min(Math.min(width, height), Math.max(1, matrixLimit)));
        List<String> sizes = new ArrayList<>();
        for (int size = 8; size <= max; size *= 2) {
            sizes.add(Integer.toString(size));
        }
        return sizes.toArray(new String[0]);
    }

    private int largestPowerOfTwoAtMost(int value) {
        int power = 1;
        while (power * 2 <= value) {
            power *= 2;
        }
        return power;
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

    private static class RoiPreviewController implements TextListener, ItemListener {
        private final ImagePlus imp;
        private final Choice roiSizeChoice;
        private final TextField matrixField;
        private final TextField centerXField;
        private final TextField centerYField;
        private final Overlay originalOverlay;
        private boolean updatingChoice = false;

        RoiPreviewController(ImagePlus imp, Choice roiSizeChoice, TextField matrixField, TextField centerXField, TextField centerYField) {
            this.imp = imp;
            this.roiSizeChoice = roiSizeChoice;
            this.matrixField = matrixField;
            this.centerXField = centerXField;
            this.centerYField = centerYField;
            this.originalOverlay = imp.getOverlay();
        }

        void install() {
            roiSizeChoice.addItemListener(this);
            matrixField.addTextListener(this);
            centerXField.addTextListener(this);
            centerYField.addTextListener(this);
            updateFromFields();
        }

        void uninstall() {
            roiSizeChoice.removeItemListener(this);
            matrixField.removeTextListener(this);
            centerXField.removeTextListener(this);
            centerYField.removeTextListener(this);
            imp.setOverlay(originalOverlay);
            imp.updateAndDraw();
        }

        @Override
        public void textValueChanged(TextEvent e) {
            updateFromFields();
        }

        @Override
        public void itemStateChanged(ItemEvent e) {
            updateFromFields();
        }

        void updateFromFields() {
            updateRoiSizeChoices();
            Rectangle rect = currentRectangle();
            if (rect == null) {
                imp.setOverlay(originalOverlay);
                imp.updateAndDraw();
                return;
            }

            boolean inside = rect.x >= 0
                && rect.y >= 0
                && rect.x + rect.width <= imp.getWidth()
                && rect.y + rect.height <= imp.getHeight();
            Roi roi = new Roi(rect.x, rect.y, rect.width, rect.height);
            roi.setName("NPS ROI preview");
            roi.setStrokeColor(inside ? Color.YELLOW : Color.RED);
            roi.setStrokeWidth(2.0);
            Overlay overlay = originalOverlay == null ? new Overlay() : originalOverlay.duplicate();
            addThumbnailPreview(overlay, rect, inside);
            overlay.add(roi);
            imp.setOverlay(overlay);
            imp.updateAndDraw();
        }

        private void updateRoiSizeChoices() {
            if (updatingChoice) {
                return;
            }
            int matrixLimit = currentMatrixLimit();
            if (matrixLimit < 8) {
                return;
            }
            String[] allowed = allowedPowerOfTwoRoiSizes(imp.getWidth(), imp.getHeight(), matrixLimit);
            if (allowed.length == 0) {
                return;
            }
            String selected = roiSizeChoice.getSelectedItem();
            if (choiceItemsEqual(roiSizeChoice, allowed)) {
                return;
            }
            updatingChoice = true;
            roiSizeChoice.removeAll();
            int selectedValue = parsePositiveInt(selected, Integer.parseInt(allowed[allowed.length - 1]));
            String replacement = allowed[0];
            for (String value : allowed) {
                roiSizeChoice.add(value);
                if (Integer.parseInt(value) <= selectedValue) {
                    replacement = value;
                }
            }
            roiSizeChoice.select(replacement);
            updatingChoice = false;
        }

        private int currentMatrixLimit() {
            return parsePositiveInt(matrixField.getText(), Math.min(imp.getWidth(), imp.getHeight()));
        }

        private boolean choiceItemsEqual(Choice choice, String[] allowed) {
            if (choice.getItemCount() != allowed.length) {
                return false;
            }
            for (int i = 0; i < allowed.length; i++) {
                if (!allowed[i].equals(choice.getItem(i))) {
                    return false;
                }
            }
            return true;
        }

        private int parsePositiveInt(String text, int fallback) {
            try {
                int value = (int) Math.round(Double.parseDouble(text));
                return value > 0 ? value : fallback;
            } catch (NumberFormatException e) {
                return fallback;
            }
        }

        private static String[] allowedPowerOfTwoRoiSizes(int width, int height, int matrixLimit) {
            int max = 1;
            int bound = Math.min(Math.min(width, height), Math.max(1, matrixLimit));
            while (max * 2 <= bound) {
                max *= 2;
            }
            List<String> sizes = new ArrayList<>();
            for (int size = 8; size <= max; size *= 2) {
                sizes.add(Integer.toString(size));
            }
            return sizes.toArray(new String[0]);
        }

        private void addThumbnailPreview(Overlay overlay, Rectangle rect, boolean inside) {
            Rectangle imageBounds = new Rectangle(0, 0, imp.getWidth(), imp.getHeight());
            Rectangle cropBounds = rect.intersection(imageBounds);
            if (cropBounds.isEmpty()) {
                return;
            }

            int thumbnailSize = Math.max(32, Math.min(128, Math.min(imp.getWidth(), imp.getHeight()) / 5));
            Rectangle thumbnailBounds = thumbnailBounds(thumbnailSize, rect);
            try {
                ImageProcessor processor = imp.getProcessor().duplicate();
                processor.setMinAndMax(imp.getDisplayRangeMin(), imp.getDisplayRangeMax());
                processor.setRoi(cropBounds);
                ImageProcessor thumbnail = processor.crop().resize(thumbnailSize, thumbnailSize, true).convertToByte(true).convertToRGB();
                ImageRoi thumbnailRoi = new ImageRoi(thumbnailBounds.x, thumbnailBounds.y, thumbnail);
                thumbnailRoi.setName("NPS ROI thumbnail");
                overlay.add(thumbnailRoi);

                Roi border = new Roi(thumbnailBounds.x, thumbnailBounds.y, thumbnailBounds.width, thumbnailBounds.height);
                border.setName("NPS ROI thumbnail border");
                border.setStrokeColor(inside ? Color.YELLOW : Color.RED);
                border.setStrokeWidth(2.0);
                overlay.add(border);
            } catch (RuntimeException e) {
                // Keep the ROI rectangle preview even if a display-specific thumbnail conversion fails.
            }
        }

        private Rectangle thumbnailBounds(int thumbnailSize, Rectangle selectedRoi) {
            int margin = Math.max(4, thumbnailSize / 12);
            Rectangle[] candidates = new Rectangle[] {
                new Rectangle(imp.getWidth() - thumbnailSize - margin, margin, thumbnailSize, thumbnailSize),
                new Rectangle(margin, margin, thumbnailSize, thumbnailSize),
                new Rectangle(imp.getWidth() - thumbnailSize - margin, imp.getHeight() - thumbnailSize - margin, thumbnailSize, thumbnailSize),
                new Rectangle(margin, imp.getHeight() - thumbnailSize - margin, thumbnailSize, thumbnailSize)
            };
            for (Rectangle candidate : candidates) {
                if (!candidate.intersects(selectedRoi)) {
                    return candidate;
                }
            }
            return candidates[0];
        }

        private Rectangle currentRectangle() {
            try {
                int size = Integer.parseInt(roiSizeChoice.getSelectedItem());
                double centerX = Double.parseDouble(centerXField.getText());
                double centerY = Double.parseDouble(centerYField.getText());
                int x0 = (int) Math.round(centerX - size / 2.0);
                int y0 = (int) Math.round(centerY - size / 2.0);
                return new Rectangle(x0, y0, size, size);
            } catch (NumberFormatException e) {
                return null;
            }
        }
    }

    private static class ImageCoordinateProbe extends MouseAdapter {
        private final ImagePlus imp;
        private final ImageCanvas canvas;
        private final Calibration cal;
        private final TextField centerXField;
        private final TextField centerYField;
        private final RoiPreviewController preview;

        ImageCoordinateProbe(ImagePlus imp, ImageCanvas canvas, Calibration cal, TextField centerXField, TextField centerYField, RoiPreviewController preview) {
            this.imp = imp;
            this.canvas = canvas;
            this.cal = cal;
            this.centerXField = centerXField;
            this.centerYField = centerYField;
            this.preview = preview;
        }

        void uninstall() {
            canvas.removeMouseListener(this);
            canvas.removeMouseMotionListener(this);
        }

        @Override
        public void mouseMoved(MouseEvent e) {
            report(e, false);
        }

        @Override
        public void mouseClicked(MouseEvent e) {
            if (e.isPopupTrigger() || e.getButton() != MouseEvent.BUTTON1) {
                return;
            }
            report(e, true);
        }

        private void report(MouseEvent e, boolean setCenter) {
            int x = canvas.offScreenX(e.getX());
            int y = canvas.offScreenY(e.getY());
            if (x < 0 || y < 0 || x >= imp.getWidth() || y >= imp.getHeight()) {
                return;
            }
            double raw = imp.getProcessor().getf(x, y);
            double value = cal != null ? cal.getCValue(raw) : raw;
            if (setCenter) {
                centerXField.setText(Integer.toString(x));
                centerYField.setText(Integer.toString(y));
                if (preview != null) {
                    preview.updateFromFields();
                }
            }
            String message = String.format(
                Locale.US,
                "NPS ROI probe: x=%d, y=%d, value=%.3f%s",
                x,
                y,
                value,
                setCenter ? " -> ROI center" : ""
            );
            IJ.showStatus(message);
            if (setCenter) {
                IJ.log(message);
            }
        }
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

    private static class NpsParameters {
        final double pixelSizeMm;
        final double fovMm;
        final int matrix;
        final boolean geometryDicomDerived;
        final String geometrySource;
        final Rectangle roiBounds;
        final String trendRemoval;
        final boolean useHannWindow;

        NpsParameters(double pixelSizeMm, double fovMm, int matrix, boolean geometryDicomDerived,
                      String geometrySource, Rectangle roiBounds, String trendRemoval, boolean useHannWindow) {
            this.pixelSizeMm = pixelSizeMm;
            this.fovMm = fovMm;
            this.matrix = matrix;
            this.geometryDicomDerived = geometryDicomDerived;
            this.geometrySource = geometrySource;
            this.roiBounds = roiBounds;
            this.trendRemoval = trendRemoval;
            this.useHannWindow = useHannWindow;
        }
    }

    private static class NpsResult {
        final double[] frequencies;
        final double[][] perSlice;
        final double[] average;
        final double frequencyStep;
        final double nyquist;
        final NpsParameters parameters;

        NpsResult(double[] frequencies, double[][] perSlice, double[] average,
                  double frequencyStep, double nyquist, NpsParameters parameters) {
            this.frequencies = frequencies;
            this.perSlice = perSlice;
            this.average = average;
            this.frequencyStep = frequencyStep;
            this.nyquist = nyquist;
            this.parameters = parameters;
        }
    }
}
