import ij.IJ;
import ij.ImagePlus;
import ij.WindowManager;
import ij.gui.GenericDialog;
import ij.gui.ImageCanvas;
import ij.gui.ImageRoi;
import ij.gui.NonBlockingGenericDialog;
import ij.gui.Overlay;
import ij.gui.OvalRoi;
import ij.gui.Plot;
import ij.gui.PlotWindow;
import ij.gui.Roi;
import ij.gui.YesNoCancelDialog;
import ij.io.SaveDialog;
import ij.measure.Calibration;
import ij.measure.CurveFitter;
import ij.measure.ResultsTable;
import ij.io.FileInfo;
import ij.plugin.PlugIn;
import ij.process.ColorProcessor;
import ij.process.FloatProcessor;
import ij.process.ImageProcessor;

import java.awt.Color;
import java.awt.Rectangle;
import java.awt.TextField;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.awt.event.TextEvent;
import java.awt.event.TextListener;
import java.util.ArrayList;
import java.util.List;
import java.util.Locale;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;

/**
 * Circular edge alignment and TTF measurement for ImageJ 1.x.
 *
 * <p>The workflow is based on the previous CircularEdgeCorrector plugin:
 * sample radial profiles around a circular insert, estimate the edge position
 * for each angle, align the profiles, average the corrected ESF, and export the
 * profile. This release adds LSF/TTF calculation from the aligned ESF.</p>
 */
public class Circular_Edge_TTF implements PlugIn {
    private static final String VERSION = "2026-05-22-fix13";
    protected static final int DEFAULT_CENTER_X = 391;
    protected static final int DEFAULT_CENTER_Y = 386;
    protected static final int DEFAULT_MAX_RADIUS = 60;
    protected static final int DEFAULT_ANGLE_SAMPLES = 360;
    protected static final double OUTER_MIN_HU = -50.0;
    protected static final double DEFAULT_DIAMETER_MM = 28.0;
    protected static final double DEFAULT_TOLERANCE_MM = 5.0;
    protected static final double DEFAULT_PIXEL_SIZE_MM = 1.0;
    private static final Pattern NUMBER_PATTERN = Pattern.compile("[-+]?\\d*\\.?\\d+(?:[Ee][-+]?\\d+)?");

    protected double centerX = DEFAULT_CENTER_X;
    protected double centerY = DEFAULT_CENTER_Y;
    protected double objectDiameterMm = DEFAULT_DIAMETER_MM;
    protected double toleranceMm = DEFAULT_TOLERANCE_MM;
    protected double pixelSizeMm = DEFAULT_PIXEL_SIZE_MM;
    protected double dicomFovMm = Double.NaN;
    protected int dicomMatrix = 0;
    protected boolean geometryDicomDerived = false;
    protected String pixelSizeSource = "not calculated";
    protected int maxRadius = DEFAULT_MAX_RADIUS;
    protected int angleSamples = DEFAULT_ANGLE_SAMPLES;
    protected int fftSamples = 64;
    protected boolean exportCsv = true;

    protected double expectedEdgePx = -1.0;
    protected double tolerancePx = 0.0;

    @Override
    public void run(String arg) {
        IJ.log("Circular Edge TTF " + VERSION + " loaded as class " + getClass().getName());
        ImagePlus imp = IJ.getImage();
        if (imp == null) {
            IJ.error("No image", "Open an image first.");
            return;
        }

        Calibration cal = imp.getCalibration();
        int width = imp.getWidth();
        int height = imp.getHeight();
        PixelGeometry geometry = inferPixelGeometry(imp, cal, width, height);
        if (!geometry.isDicomDerived) {
            IJ.log("Circular Edge TTF: DICOM FOV/pixel-spacing tags were not found. Source attempted: " + geometry.source);
            geometry = requestManualGeometry(geometry, width);
            if (geometry == null) {
                IJ.log("Circular Edge TTF cancelled: DICOM geometry missing and manual geometry input was cancelled.");
                return;
            }
        }
        pixelSizeMm = geometry.pixelSizeMm;
        dicomFovMm = geometry.fovMm;
        dicomMatrix = geometry.matrix;
        geometryDicomDerived = geometry.isDicomDerived;
        pixelSizeSource = geometry.source;
        IJ.log(String.format(
            Locale.US,
            "Circular Edge TTF geometry: pixelSize=%.9f mm/pixel, FOV=%.6f mm, matrix=%d, dicomDerived=%s, source=%s",
            pixelSizeMm,
            dicomFovMm,
            dicomMatrix,
            geometryDicomDerived,
            pixelSizeSource
        ));

        if (!readParameters(imp, width, height)) {
            return;
        }

        expectedEdgePx = (objectDiameterMm / 2.0) / pixelSizeMm;
        tolerancePx = Math.max(1.0, toleranceMm / pixelSizeMm);
        int radialLength = maxRadius + 1;
        int startSlice = Math.max(1, imp.getCurrentSlice());
        int endSlice = imp.getStackSize();

        PlotWindow[] ttfWindows = null;
        for (int slice = startSlice; slice <= endSlice; slice++) {
            ImageProcessor ip = imp.getStack().getProcessor(slice).duplicate();
            if (imp.getType() == ImagePlus.COLOR_RGB || ip instanceof ColorProcessor) {
                ip = ip.convertToByte(true);
            }
            ip.setInterpolationMethod(ImageProcessor.BILINEAR);

            SliceResult autoResult = analyzeSlice(ip, cal, width, height, radialLength);
            autoResult.rawPolar.resetMinAndMax();
            autoResult.alignedPolar.resetMinAndMax();

            showOrReplaceImage("Polar Raw (HU)", autoResult.rawPolar);
            showOrReplaceImage("Aligned Polar (HU)", autoResult.alignedPolar);

            double[] radiusPx = new double[radialLength];
            double[] radiusMm = new double[radialLength];
            for (int r = 0; r < radialLength; r++) {
                radiusPx[r] = r;
                radiusMm[r] = r * pixelSizeMm;
            }

            SliceResult tuned = adjustPlateausInteractively(autoResult, radiusPx);
            if (tuned == null) {
                IJ.log(String.format("Processing cancelled at slice %d", slice));
                break;
            }

            double[] segmentedFull = buildSegmentedProfile(
                tuned.meanProfile,
                tuned.highPlateau,
                tuned.lowPlateau,
                tuned.highIntersection,
                tuned.lowIntersection
            );
            TtfResult ttf = computeTtf(segmentedFull, pixelSizeMm);
            ttfWindows = showTtfPlots(ttf, ttfWindows);

            boolean csvSaved = false;
            if (exportCsv) {
                csvSaved = exportAnalysisToCSV(slice, radiusPx, radiusMm, tuned.meanProfile, segmentedFull, ttf);
            }

            ResultsTable rt = ResultsTable.getResultsTable();
            if (rt == null) {
                rt = new ResultsTable();
            }
            rt.incrementCounter();
            rt.addValue("Slice", slice);
            rt.addValue("HighPlateau_HU", tuned.highPlateau);
            rt.addValue("LowPlateau_HU", tuned.lowPlateau);
            rt.addValue("HighIntersection_px", tuned.highIntersection);
            rt.addValue("LowIntersection_px", tuned.lowIntersection);
            rt.addValue("PixelSize_mm", pixelSizeMm);
            rt.addValue("Geometry_FOV_mm", dicomFovMm);
            rt.addValue("Geometry_Matrix", dicomMatrix);
            rt.addValue("DICOM_FOV_mm", dicomFovMm);
            rt.addValue("DICOM_Matrix", dicomMatrix);
            rt.addValue("Geometry_DICOM_Derived", geometryDicomDerived ? 1 : 0);
            rt.addValue("TTF_FFT_N", ttf.fftLength);
            rt.addValue("TTF_FrequencyStep_cycles_per_mm", ttf.frequencyStep);
            rt.addValue("TTF50_cycles_per_mm", ttf.ttf50);
            rt.addValue("TTF10_cycles_per_mm", ttf.ttf10);
            rt.show("Results");

            IJ.log(String.format(
                "Slice %d: HighPlateau=%.2f HU LowPlateau=%.2f HU HighInt=%.2f px LowInt=%.2f px pixel=%.9f mm FOV=%.6f mm matrix=%d FFT_N=%d df=%.6f cy/mm TTF50=%.4f cy/mm TTF10=%.4f cy/mm",
                slice,
                tuned.highPlateau,
                tuned.lowPlateau,
                tuned.highIntersection,
                tuned.lowIntersection,
                pixelSizeMm,
                dicomFovMm,
                dicomMatrix,
                ttf.fftLength,
                ttf.frequencyStep,
                ttf.ttf50,
                ttf.ttf10
            ));

            if (confirmFinishAfterCalculation("Circular Edge TTF", csvSaved)) {
                break;
            }
        }

        closePlots(ttfWindows);
        closeWindowByTitle("Flattened Profile");
        closeWindowByTitle("Segmented Profile");
        closeWindowByTitle("Polar Raw (HU)");
        closeWindowByTitle("Aligned Polar (HU)");
    }

    private boolean readParameters(ImagePlus imp, int width, int height) {
        GenericDialog gd = new NonBlockingGenericDialog("Circular Edge TTF Parameters");
        gd.addMessage(String.format(
            Locale.US,
            "%s geometry: FOV %.6f mm / matrix %d = %.9f mm/pixel%nSource: %s",
            geometryDicomDerived ? "DICOM-derived" : "Manual",
            dicomFovMm,
            dicomMatrix,
            pixelSizeMm,
            pixelSizeSource
        ));
        gd.addNumericField("Center X (pixels)", centerX, 2);
        gd.addNumericField("Center Y (pixels)", centerY, 2);
        gd.addNumericField("Object diameter (mm)", objectDiameterMm, 2);
        gd.addNumericField("Edge tolerance (mm)", toleranceMm, 2);
        gd.addNumericField("Sampling radius (pixels)", maxRadius, 0);
        gd.addNumericField("Angular samples", angleSamples, 0);
        gd.addNumericField("TTF FFT samples", fftSamples, 0);
        gd.addMessage("Move the cursor over the image to read x/y/value. Click the image to set Center X/Y.");
        gd.addCheckbox("Export CSV per slice", exportCsv);
        TtfCenterPreviewController preview = installTtfCenterPreview(imp, gd);
        TtfImageCoordinateProbe probe = installTtfImageCoordinateProbe(imp, gd, preview);
        gd.showDialog();
        if (probe != null) {
            probe.uninstall();
        }
        if (preview != null) {
            preview.uninstall();
        }
        if (gd.wasCanceled()) {
            return false;
        }

        double enteredX = gd.getNextNumber();
        double enteredY = gd.getNextNumber();
        double enteredDiameter = gd.getNextNumber();
        double enteredTolerance = gd.getNextNumber();
        double enteredRadius = gd.getNextNumber();
        double enteredAngles = gd.getNextNumber();
        double enteredFftSamples = gd.getNextNumber();
        boolean enteredExport = gd.getNextBoolean();

        if (gd.invalidNumber() || Double.isNaN(enteredX) || Double.isNaN(enteredY)
            || Double.isNaN(enteredDiameter)
            || Double.isNaN(enteredTolerance) || Double.isNaN(enteredRadius)
            || Double.isNaN(enteredAngles) || Double.isNaN(enteredFftSamples)) {
            IJ.error("Invalid input", "All numeric entries must be valid numbers.");
            return false;
        }
        if (enteredX < 0 || enteredX >= width || enteredY < 0 || enteredY >= height) {
            IJ.error("Invalid center", "Center must be inside the image bounds.");
            return false;
        }
        if (enteredDiameter <= 0 || enteredTolerance <= 0) {
            IJ.error("Invalid edge geometry", "Object diameter and edge tolerance must be greater than zero.");
            return false;
        }
        if (enteredRadius < 8 || enteredAngles < 16) {
            IJ.error("Invalid sampling", "Use a sampling radius >= 8 px and angular samples >= 16.");
            return false;
        }
        if (enteredFftSamples < 8) {
            IJ.error("Invalid FFT samples", "Use at least 8 FFT samples. The validated default is 64.");
            return false;
        }

        centerX = enteredX;
        centerY = enteredY;
        objectDiameterMm = enteredDiameter;
        toleranceMm = enteredTolerance;
        maxRadius = Math.max(8, (int) Math.round(enteredRadius));
        angleSamples = Math.max(16, (int) Math.round(enteredAngles));
        fftSamples = Math.max(8, (int) Math.round(enteredFftSamples));
        exportCsv = enteredExport;
        return true;
    }

    private TtfCenterPreviewController installTtfCenterPreview(ImagePlus imp, GenericDialog gd) {
        List<?> fields = gd.getNumericFields();
        if (fields == null || fields.size() < 5) {
            return null;
        }
        TtfCenterPreviewController preview = new TtfCenterPreviewController(
            imp,
            (TextField) fields.get(0),
            (TextField) fields.get(1),
            (TextField) fields.get(2),
            (TextField) fields.get(4)
        );
        preview.install();
        return preview;
    }

    private TtfImageCoordinateProbe installTtfImageCoordinateProbe(ImagePlus imp, GenericDialog gd, TtfCenterPreviewController preview) {
        ImageCanvas canvas = imp.getCanvas();
        if (canvas == null) {
            return null;
        }
        List<?> fields = gd.getNumericFields();
        if (fields == null || fields.size() < 2) {
            return null;
        }
        TtfImageCoordinateProbe probe = new TtfImageCoordinateProbe(
            imp,
            canvas,
            imp.getCalibration(),
            (TextField) fields.get(0),
            (TextField) fields.get(1),
            preview
        );
        canvas.addMouseListener(probe);
        canvas.addMouseMotionListener(probe);
        return probe;
    }

    private double inferPixelSizeMm(Calibration cal) {
        if (cal == null) {
            return DEFAULT_PIXEL_SIZE_MM;
        }
        double candidate = cal.pixelWidth > 0 ? cal.pixelWidth : cal.pixelHeight;
        if (candidate <= 0) {
            return DEFAULT_PIXEL_SIZE_MM;
        }
        String unit = cal.getUnit();
        if (unit == null || unit.equalsIgnoreCase("pixel") || unit.equalsIgnoreCase("pixels")) {
            return DEFAULT_PIXEL_SIZE_MM;
        }
        if (unit.equalsIgnoreCase("cm")) {
            return candidate * 10.0;
        }
        if (unit.equalsIgnoreCase("um") || unit.equalsIgnoreCase("micron") || unit.equalsIgnoreCase("microns")) {
            return candidate / 1000.0;
        }
        return candidate;
    }

    protected PixelGeometry inferPixelGeometry(ImagePlus imp, Calibration cal, int width, int height) {
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

        if (rows.intValue() != height || columns.intValue() != width) {
            IJ.log(String.format(
                Locale.US,
                "Circular Edge TTF warning: Image dimensions (%dx%d) differ from DICOM Rows/Columns (%dx%d). Using DICOM Columns for frequency scale.",
                width,
                height,
                columns,
                rows
            ));
        } else if ("ImageJ image dimensions".equals(matrixSource)) {
            IJ.log("Circular Edge TTF: DICOM Rows/Columns tags were not found; using ImageJ image width/height as matrix.");
        }

        Double reconstructionDiameter = parseDicomDouble(info, "0018,1100", "Reconstruction Diameter");
        if (isPositive(reconstructionDiameter)) {
            double fovMm = reconstructionDiameter.doubleValue();
            return new PixelGeometry(
                fovMm / columns.doubleValue(),
                fovMm,
                columns.intValue(),
                true,
                String.format(Locale.US, "DICOM Reconstruction Diameter (0018,1100)=%.9g mm / matrix=%d (%s)", fovMm, columns, matrixSource)
            );
        }

        double[] pixelSpacing = parseDicomDoubles(info, "0028,0030", "Pixel Spacing");
        if (pixelSpacing.length > 0 && isPositive(pixelSpacing[0])) {
            double rowSpacing = pixelSpacing[0];
            double columnSpacing = pixelSpacing.length > 1 && isPositive(pixelSpacing[1]) ? pixelSpacing[1] : rowSpacing;
            if (!nearlyEqual(rowSpacing, columnSpacing)) {
                return calibrationFallback(
                    cal,
                    width,
                    String.format(Locale.US, "anisotropic DICOM Pixel Spacing (0028,0030)=%.9g/%.9g mm", rowSpacing, columnSpacing)
                );
            }
            double fovMm = columnSpacing * columns.doubleValue();
            return new PixelGeometry(
                fovMm / columns.doubleValue(),
                fovMm,
                columns.intValue(),
                true,
                String.format(Locale.US, "DICOM Pixel Spacing (0028,0030)=%.9g mm * matrix=%d (%s)", columnSpacing, columns, matrixSource)
            );
        }

        return calibrationFallback(cal, width, "missing DICOM Reconstruction Diameter (0018,1100) and Pixel Spacing (0028,0030)");
    }

    private PixelGeometry calibrationFallback(Calibration cal, int width, String reason) {
        double fallbackPixelSize = inferPixelSizeMm(cal);
        return new PixelGeometry(
            fallbackPixelSize,
            fallbackPixelSize * width,
            width,
            false,
            reason + "; calibration fallback would be " + String.format(Locale.US, "%.9g mm/pixel", fallbackPixelSize)
        );
    }

    protected PixelGeometry requestManualGeometry(PixelGeometry fallback, int width) {
        double defaultFov = fallback != null && isPositive(Double.valueOf(fallback.fovMm))
            ? fallback.fovMm
            : DEFAULT_PIXEL_SIZE_MM * width;
        int defaultMatrix = fallback != null && fallback.matrix > 0 ? fallback.matrix : width;

        GenericDialog gd = new GenericDialog("Manual TTF Frequency Geometry");
        gd.addMessage(
            "DICOM FOV/pixel-spacing tags were not found.\n"
            + "Enter the scan FOV and matrix size used for this image.\n"
            + "The spatial-frequency step will be calculated as matrix / (FFT_N * FOV)."
        );
        if (fallback != null && fallback.source != null) {
            gd.addMessage("Reason: " + fallback.source);
        }
        gd.addNumericField("FOV (mm)", defaultFov, 6);
        gd.addNumericField("Matrix size (pixels)", defaultMatrix, 0);
        gd.showDialog();
        if (gd.wasCanceled()) {
            return null;
        }
        double enteredFov = gd.getNextNumber();
        double enteredMatrix = gd.getNextNumber();
        if (gd.invalidNumber() || enteredFov <= 0 || enteredMatrix <= 0
            || Double.isNaN(enteredFov) || Double.isNaN(enteredMatrix)) {
            IJ.error("Invalid manual geometry", "FOV and matrix size must be greater than zero.");
            return requestManualGeometry(fallback, width);
        }

        int matrix = Math.max(1, (int) Math.round(enteredMatrix));
        return new PixelGeometry(
            enteredFov / matrix,
            enteredFov,
            matrix,
            false,
            String.format(Locale.US, "Manual FOV/matrix input: FOV=%.9g mm, matrix=%d", enteredFov, matrix)
        );
    }

    private String dicomInfo(ImagePlus imp) {
        if (imp == null) {
            return "";
        }
        StringBuilder sb = new StringBuilder();
        String info = imp.getInfoProperty();
        if (info != null && info.trim().length() > 0) {
            appendDicomInfo(sb, info);
        }
        Object property = imp.getProperty("Info");
        if (property != null) {
            appendDicomInfo(sb, property.toString());
        }
        FileInfo fileInfo = imp.getOriginalFileInfo();
        if (fileInfo != null && fileInfo.info != null) {
            appendDicomInfo(sb, fileInfo.info);
        }
        if (imp.getStack() != null) {
            int current = Math.max(1, Math.min(imp.getCurrentSlice(), imp.getStackSize()));
            appendDicomInfo(sb, imp.getStack().getSliceLabel(current));
            if (current != 1 && imp.getStackSize() >= 1) {
                appendDicomInfo(sb, imp.getStack().getSliceLabel(1));
            }
        }
        return sb.toString();
    }

    private void appendDicomInfo(StringBuilder sb, String text) {
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

        if (label == null || label.length() == 0) {
            return null;
        }
        String normalizedLabel = label.toLowerCase(Locale.US);
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

    private boolean isPositive(Number value) {
        return value != null && value.doubleValue() > 0.0 && !Double.isNaN(value.doubleValue()) && !Double.isInfinite(value.doubleValue());
    }

    private boolean nearlyEqual(double a, double b) {
        double tolerance = Math.max(1e-6, Math.max(Math.abs(a), Math.abs(b)) * 1e-4);
        return Math.abs(a - b) <= tolerance;
    }

    protected SliceResult adjustPlateausInteractively(SliceResult initial, double[] axis) {
        SliceResult current = initial;
        PlotWindow[] windows = null;
        boolean cancelled = false;
        while (true) {
            windows = showProfilePlots(current, axis, windows);

            NonBlockingGenericDialog gd = new NonBlockingGenericDialog("Adjust Plateau / Intersections");
            String message = "Inspect the displayed plots. Modify the values below if needed.";
            message += " Use Redraw to apply typed values without measuring, or Measure TTF to continue.";
            message += " If a plateau is changed and its intersection is left unchanged, the intersection is recalculated.";
            message += String.format(
                "%nFlattened Profile: click at a radius position to insert the curve CT value and intersection at that X-axis position. Inner side sets High; outer side sets Low. Shift-click=High, Ctrl/Alt-click=Low."
            );
            if (expectedEdgePx > 0) {
                message += String.format("%nExpected edge approx %.2f px (+/- %.2f)", expectedEdgePx, tolerancePx);
            }
            gd.addMessage(message);
            gd.addNumericField("High Plateau (HU)", current.highPlateau, 2);
            gd.addNumericField("Low Plateau (HU)", current.lowPlateau, 2);
            gd.addNumericField("High Intersection (radius px)", current.highIntersection, 2);
            gd.addNumericField("Low Intersection (radius px)", current.lowIntersection, 2);
            ProfileAdjustmentController controller = installProfileAdjustmentController(windows, current, axis, gd);
            if (controller != null) {
                gd.addButton("Redraw", e -> controller.redrawFromFields());
            }
            beforePlateauAdjustmentDialogShown(gd, current, axis);
            gd.setOKLabel("Measure TTF");
            gd.showDialog();
            if (controller != null) {
                windows = controller.getWindows();
                current = controller.getCurrent();
            }
            if (gd.wasCanceled() && !acceptPlateauAdjustmentDialogDismissal()) {
                cancelled = true;
                break;
            }
            double highField = gd.getNextNumber();
            double lowField = gd.getNextNumber();
            double highIntField = gd.getNextNumber();
            double lowIntField = gd.getNextNumber();

            double highPlateau = Double.isNaN(highField) ? current.highPlateau : highField;
            double lowPlateau = Double.isNaN(lowField) ? current.lowPlateau : lowField;
            double highIntersection = Double.isNaN(highIntField) ? current.highIntersection : highIntField;
            double lowIntersection = Double.isNaN(lowIntField) ? current.lowIntersection : lowIntField;

            boolean highPlateauChanged = !Double.isNaN(highField) && Math.abs(highPlateau - current.highPlateau) > 1e-9;
            boolean lowPlateauChanged = !Double.isNaN(lowField) && Math.abs(lowPlateau - current.lowPlateau) > 1e-9;
            boolean highIntersectionChanged = !Double.isNaN(highIntField) && Math.abs(highIntersection - current.highIntersection) > 1e-9;
            boolean lowIntersectionChanged = !Double.isNaN(lowIntField) && Math.abs(lowIntersection - current.lowIntersection) > 1e-9;

            if ((highPlateauChanged && !highIntersectionChanged) || (lowPlateauChanged && !lowIntersectionChanged)) {
                IntersectionPair auto = recomputeIntersections(
                    current.meanProfile,
                    highPlateau,
                    lowPlateau,
                    current.highIntersection,
                    current.lowIntersection
                );
                if (highPlateauChanged && !highIntersectionChanged) {
                    highIntersection = auto.high;
                }
                if (lowPlateauChanged && !lowIntersectionChanged) {
                    lowIntersection = auto.low;
                }
            }

            highIntersection = Math.max(0, Math.min(maxRadius, highIntersection));
            lowIntersection = Math.max(0, Math.min(maxRadius, lowIntersection));
            if (lowIntersection < highIntersection) {
                double tmp = highIntersection;
                highIntersection = lowIntersection;
                lowIntersection = tmp;
            }

            current = new SliceResult(
                current.rawPolar,
                current.alignedPolar,
                current.meanProfile,
                highPlateau,
                lowPlateau,
                highIntersection,
                lowIntersection
            );
            break;
        }
        closePlots(windows);
        return cancelled ? null : current;
    }

    protected void beforePlateauAdjustmentDialogShown(GenericDialog gd, SliceResult current, double[] axis) {
        // Subclasses can add workflow-specific buttons without changing the single-slice plugin.
    }

    protected boolean acceptPlateauAdjustmentDialogDismissal() {
        return false;
    }

    private PlotWindow[] showProfilePlots(SliceResult result, double[] axis, PlotWindow[] previous) {
        closePlots(previous);

        EdgeSegmentData segData = extractSegmentData(result);
        PlotLimits sharedLimits = calculateProfilePlotLimits(result, segData);

        Plot flattened = new Plot("Flattened Profile", "Radius (px)", "Intensity (HU)");
        flattened.setFrameSize(720, 360);
        flattened.setMaxIntervals(14);
        flattened.setLimits(sharedLimits.xMin, sharedLimits.xMax, sharedLimits.yMin, sharedLimits.yMax);
        flattened.setColor(Color.RED);
        flattened.add("line", axis, result.meanProfile);
        flattened.setColor(Color.BLUE);
        flattened.add("line", new double[]{0, maxRadius}, new double[]{result.highPlateau, result.highPlateau});
        flattened.add("line", new double[]{0, maxRadius}, new double[]{result.lowPlateau, result.lowPlateau});
        flattened.setColor(Color.WHITE);
        flattened.addPoints(new double[]{result.highIntersection}, new double[]{result.highPlateau}, Plot.CIRCLE);
        flattened.addPoints(new double[]{result.lowIntersection}, new double[]{result.lowPlateau}, Plot.CIRCLE);
        PlotWindow flatWin = flattened.show();
        if (flatWin != null) {
            flatWin.setLocation(100, 120);
            flatWin.toFront();
        }

        Plot segmented = new Plot("Segmented Profile", "Radius (px)", "Intensity (HU)");
        segmented.setFrameSize(720, 360);
        segmented.setMaxIntervals(14);
        segmented.setLimits(sharedLimits.xMin, sharedLimits.xMax, sharedLimits.yMin, sharedLimits.yMax);
        segmented.setColor(Color.RED);
        segmented.setLineWidth(2);
        segmented.add("line", segData.radius, segData.intensity);
        segmented.setLineWidth(1);
        segmented.setColor(Color.BLUE);
        segmented.add("line", new double[]{0, result.highIntersection}, new double[]{result.highPlateau, result.highPlateau});
        segmented.add("line", new double[]{result.lowIntersection, maxRadius}, new double[]{result.lowPlateau, result.lowPlateau});
        segmented.setColor(Color.WHITE);
        segmented.addPoints(new double[]{result.highIntersection}, new double[]{result.highPlateau}, Plot.CIRCLE);
        segmented.addPoints(new double[]{result.lowIntersection}, new double[]{result.lowPlateau}, Plot.CIRCLE);
        PlotWindow segWin = segmented.show();
        if (segWin != null) {
            segWin.setLocation(flatWin != null ? flatWin.getX() + flatWin.getWidth() + 40 : 500, 120);
            segWin.toFront();
        }

        return new PlotWindow[]{flatWin, segWin};
    }

    private PlotLimits calculateProfilePlotLimits(SliceResult result, EdgeSegmentData segData) {
        double yMin = Math.min(result.highPlateau, result.lowPlateau);
        double yMax = Math.max(result.highPlateau, result.lowPlateau);
        int last = Math.min(maxRadius, result.meanProfile.length - 1);
        for (int r = 0; r <= last; r++) {
            double v = result.meanProfile[r];
            if (!Double.isNaN(v) && !Double.isInfinite(v)) {
                if (v < yMin) yMin = v;
                if (v > yMax) yMax = v;
            }
        }
        if (segData != null && segData.intensity != null) {
            for (double v : segData.intensity) {
                if (!Double.isNaN(v) && !Double.isInfinite(v)) {
                    if (v < yMin) yMin = v;
                    if (v > yMax) yMax = v;
                }
            }
        }
        double span = yMax - yMin;
        double pad = span < 1e-3 ? Math.max(1.0, Math.abs(yMax) * 0.05 + 1.0) : span * 0.1;
        return new PlotLimits(0.0, maxRadius, yMin - pad, yMax + pad);
    }

    private ProfileAdjustmentController installProfileAdjustmentController(PlotWindow[] windows, SliceResult current, double[] axis, GenericDialog gd) {
        if (windows == null || windows.length == 0 || windows[0] == null || current == null || gd == null) {
            return null;
        }
        java.util.Vector fields = gd.getNumericFields();
        if (fields == null || fields.size() < 4) {
            return null;
        }
        Object high = fields.get(0);
        Object low = fields.get(1);
        Object highInt = fields.get(2);
        Object lowInt = fields.get(3);
        if (!(high instanceof TextField) || !(low instanceof TextField)
            || !(highInt instanceof TextField) || !(lowInt instanceof TextField)) {
            return null;
        }
        ProfileAdjustmentController controller = new ProfileAdjustmentController(
            axis,
            current,
            windows,
            (TextField) high,
            (TextField) low,
            (TextField) highInt,
            (TextField) lowInt
        );
        controller.install();
        return controller;
    }

    protected PlotWindow[] showTtfPlots(TtfResult ttf, PlotWindow[] previous) {
        closePlots(previous);

        Plot lsfPlot = new Plot("Circular Edge LSF", "Radius (mm)", "LSF (HU/mm)");
        lsfPlot.setColor(Color.RED);
        lsfPlot.add("line", ttf.lsfRadiusMm, ttf.lsf);
        PlotWindow lsfWin = lsfPlot.show();
        if (lsfWin != null) {
            lsfWin.setLocation(100, 520);
            lsfWin.toFront();
        }

        Plot ttfPlot = new Plot("Circular Edge TTF", "Spatial frequency (cycles/mm)", "TTF");
        ttfPlot.setLimits(0, ttf.nyquist, 0, 1.1);
        ttfPlot.setColor(Color.BLUE);
        ttfPlot.add("line", ttf.frequency, ttf.ttf);
        ttfPlot.setColor(Color.GRAY);
        ttfPlot.add("line", new double[]{0, ttf.nyquist}, new double[]{0.5, 0.5});
        PlotWindow ttfWin = ttfPlot.show();
        if (ttfWin != null) {
            ttfWin.setLocation(lsfWin != null ? lsfWin.getX() + lsfWin.getWidth() + 40 : 500, 520);
            ttfWin.toFront();
        }
        return new PlotWindow[]{lsfWin, ttfWin};
    }

    protected void closePlots(PlotWindow[] windows) {
        if (windows == null) {
            return;
        }
        for (PlotWindow pw : windows) {
            if (pw != null && !pw.isClosed()) {
                pw.close();
            }
        }
    }

    protected void showOrReplaceImage(String title, FloatProcessor fp) {
        closeWindowByTitle(title);
        new ImagePlus(title, fp).show();
    }

    protected void closeWindowByTitle(String title) {
        ImagePlus existing = WindowManager.getImage(title);
        if (existing != null) {
            existing.close();
        }
    }

    protected SliceResult analyzeSlice(ImageProcessor ip, Calibration cal, int width, int height, int radialLength) {
        FloatProcessor rawPolar = new FloatProcessor(angleSamples, radialLength);
        for (int ang = 0; ang < angleSamples; ang++) {
            double theta = ang * 2.0 * Math.PI / angleSamples;
            double cosT = Math.cos(theta);
            double sinT = Math.sin(theta);
            for (int r = 0; r <= maxRadius; r++) {
                double x = centerX + r * cosT;
                double y = centerY + r * sinT;
                double value = 0;
                if (x >= 0 && y >= 0 && x < width && y < height) {
                    double raw = ip.getInterpolatedPixel(x, y);
                    value = cal != null ? cal.getCValue(raw) : raw;
                }
                rawPolar.setf(ang, r, (float) value);
            }
        }

        double initialEdge = expectedEdgePx > 0
            ? Math.max(1.0, Math.min(maxRadius - 1.0, expectedEdgePx))
            : Math.max(1.0, maxRadius * 0.75);
        float[] offsets = fitEdgeProfiles(rawPolar, initialEdge);
        double sumOff = 0.0;
        for (float off : offsets) {
            sumOff += off;
        }
        float meanOff = (float) (sumOff / offsets.length);
        for (int i = 0; i < offsets.length; i++) {
            offsets[i] -= meanOff;
        }

        FloatProcessor alignedPolar = alignPolarProfiles(rawPolar, offsets);

        double[] meanProf = new double[radialLength];
        for (int r = 0; r <= maxRadius; r++) {
            double s = 0;
            for (int ang = 0; ang < angleSamples; ang++) {
                s += alignedPolar.getf(ang, r);
            }
            meanProf[r] = s / angleSamples;
        }
        boolean droppingEdge = meanProf[0] > meanProf[maxRadius];

        int searchLower = (expectedEdgePx > 0 && tolerancePx > 0)
            ? Math.max(0, (int) Math.floor(expectedEdgePx - tolerancePx))
            : 0;
        int searchUpper = (expectedEdgePx > 0 && tolerancePx > 0)
            ? Math.min(maxRadius, (int) Math.ceil(expectedEdgePx + tolerancePx))
            : maxRadius;
        if (searchUpper <= searchLower) {
            searchLower = 0;
            searchUpper = maxRadius;
        }

        int peakIdx = searchLower;
        double maxDiff = -1;
        int gradientUpper = Math.min(maxRadius, searchUpper);
        for (int r = searchLower; r < gradientUpper; r++) {
            double d = droppingEdge ? -(meanProf[r + 1] - meanProf[r]) : (meanProf[r + 1] - meanProf[r]);
            if (d > maxDiff) {
                maxDiff = d;
                peakIdx = r;
            }
        }
        if (maxDiff < 0) {
            peakIdx = 0;
            maxDiff = -1;
            for (int r = 0; r < maxRadius; r++) {
                double d = droppingEdge ? -(meanProf[r + 1] - meanProf[r]) : (meanProf[r + 1] - meanProf[r]);
                if (d > maxDiff) {
                    maxDiff = d;
                    peakIdx = r;
                }
            }
        }

        int outerLimit = Math.min(maxRadius, Math.max(peakIdx + 1, gradientUpper));
        while (outerLimit > peakIdx && meanProf[outerLimit] < OUTER_MIN_HU) {
            outerLimit--;
        }
        if (outerLimit <= peakIdx) {
            outerLimit = Math.min(maxRadius, peakIdx + (int) Math.round(Math.max(1.0, tolerancePx)));
        }
        if (outerLimit <= peakIdx) {
            outerLimit = Math.min(maxRadius, peakIdx + 1);
        }

        int innerMid = (searchLower > 0)
            ? Math.max(0, Math.min(maxRadius, Math.round(searchLower / 2.0f)))
            : Math.max(0, Math.min(maxRadius, Math.round(peakIdx / 2.0f)));
        int outerMidBase = searchUpper > searchLower ? searchUpper : outerLimit;
        int outerMid = Math.max(0, Math.min(outerLimit, Math.round((peakIdx + outerMidBase) / 2.0f)));
        int win = 5;
        double sumH = 0;
        int cntH = 0;
        double sumL = 0;
        int cntL = 0;
        for (int r = innerMid - win; r <= innerMid + win; r++) {
            if (r >= 0 && r <= maxRadius) {
                sumH += meanProf[r];
                cntH++;
            }
        }
        for (int r = outerMid - win; r <= outerMid + win; r++) {
            if (r >= 0 && r <= outerLimit && meanProf[r] >= OUTER_MIN_HU) {
                sumL += meanProf[r];
                cntL++;
            }
        }
        double highPlateau = cntH > 0 ? sumH / cntH : meanProf[innerMid];
        double lowPlateau;
        if (cntL > 0) {
            lowPlateau = sumL / cntL;
        } else {
            int fallbackIdx = Math.max(0, Math.min(outerLimit, outerMid));
            lowPlateau = Math.max(OUTER_MIN_HU, meanProf[fallbackIdx]);
        }

        double highInt = findInnerPlateauEdgeCrossing(meanProf, highPlateau, peakIdx, 0, droppingEdge);
        double lowInt = findDirectedCrossing(meanProf, lowPlateau, peakIdx, outerLimit, outerLimit, droppingEdge);
        highInt = Math.max(0, Math.min(maxRadius, highInt));
        lowInt = Math.max(0, Math.min(maxRadius, lowInt));

        return new SliceResult(rawPolar, alignedPolar, meanProf, highPlateau, lowPlateau, highInt, lowInt);
    }

    protected double[] buildSegmentedProfile(double[] profile, double highPlateau, double lowPlateau, double highIntersection, double lowIntersection) {
        double[] segmented = new double[profile.length];
        for (int r = 0; r < segmented.length; r++) {
            if (r <= highIntersection) {
                segmented[r] = highPlateau;
            } else if (r >= lowIntersection) {
                segmented[r] = lowPlateau;
            } else {
                segmented[r] = profile[r];
            }
        }
        return segmented;
    }

    private IntersectionPair recomputeIntersections(double[] profile, double highPlateau, double lowPlateau,
                                                    double fallbackHigh, double fallbackLow) {
        if (profile == null || profile.length < 2) {
            return new IntersectionPair(fallbackHigh, fallbackLow);
        }
        int last = Math.min(maxRadius, profile.length - 1);
        int peakIdx = findEdgePeakIndex(profile, last);
        boolean droppingEdge = profile[0] > profile[last];

        double high = findInnerPlateauEdgeCrossing(profile, highPlateau, peakIdx, fallbackHigh, droppingEdge);
        double low = findDirectedCrossing(profile, lowPlateau, peakIdx, last, fallbackLow, droppingEdge);
        return new IntersectionPair(high, low);
    }

    private int findEdgePeakIndex(double[] profile, int last) {
        int searchLower = (expectedEdgePx > 0 && tolerancePx > 0)
            ? Math.max(0, (int) Math.floor(expectedEdgePx - tolerancePx))
            : 0;
        int searchUpper = (expectedEdgePx > 0 && tolerancePx > 0)
            ? Math.min(last, (int) Math.ceil(expectedEdgePx + tolerancePx))
            : last;
        if (searchUpper <= searchLower) {
            searchLower = 0;
            searchUpper = last;
        }

        boolean droppingEdge = profile[0] > profile[last];
        int bestIdx = searchLower;
        double bestDiff = -1.0;
        for (int r = searchLower; r < searchUpper; r++) {
            double diff = droppingEdge ? -(profile[r + 1] - profile[r]) : (profile[r + 1] - profile[r]);
            if (diff > bestDiff) {
                bestDiff = diff;
                bestIdx = r;
            }
        }
        return bestIdx;
    }

    private double findDirectedCrossing(double[] profile, double target, int start, int end,
                                        double fallback, boolean droppingEdge) {
        if (Double.isNaN(target) || end <= start) {
            return fallback;
        }
        for (int r = start; r < end; r++) {
            if (isDirectedCrossing(profile[r], profile[r + 1], target, droppingEdge)) {
                return crossingRadius(r, profile[r], profile[r + 1], target);
            }
        }
        return fallback;
    }

    private double findInnerPlateauEdgeCrossing(double[] profile, double target, int edgeIdx,
                                                double fallback, boolean droppingEdge) {
        if (profile == null || profile.length < 2 || Double.isNaN(target) || edgeIdx <= 0) {
            return fallback;
        }
        int start = Math.min(edgeIdx - 1, profile.length - 2);
        for (int r = start; r >= 0; r--) {
            if (isDirectedCrossing(profile[r], profile[r + 1], target, droppingEdge)) {
                return crossingRadius(r, profile[r], profile[r + 1], target);
            }
        }
        return fallback;
    }

    private boolean isDirectedCrossing(double v0, double v1, double target, boolean droppingEdge) {
        return droppingEdge
            ? (v0 >= target && v1 <= target)
            : (v0 <= target && v1 >= target);
    }

    private double crossingRadius(int r, double v0, double v1, double target) {
        double denom = v1 - v0;
        if (Math.abs(denom) < 1e-12) {
            return r;
        }
        double frac = (target - v0) / denom;
        frac = Math.max(0.0, Math.min(1.0, frac));
        return r + frac;
    }

    protected TtfResult computeTtf(double[] esf, double spacingMm) {
        int fftLength = Math.max(8, fftSamples);
        double[] lsf = buildPointSpread64Style(esf, spacingMm, fftLength);
        double[] lsfRadiusMm = new double[fftLength];
        for (int i = 0; i < fftLength; i++) {
            lsfRadiusMm[i] = i * spacingMm;
        }

        int bins = fftLength / 2 + 1;
        double[] frequency = new double[bins];
        double[] ttf = new double[bins];
        double zeroMagnitude = 0.0;
        for (int k = 0; k < bins; k++) {
            double re = 0.0;
            double im = 0.0;
            for (int n = 0; n < fftLength; n++) {
                double angle = -2.0 * Math.PI * k * n / fftLength;
                re += lsf[n] * Math.cos(angle);
                im += lsf[n] * Math.sin(angle);
            }
            double mag = Math.sqrt(re * re + im * im);
            if (k == 0) {
                zeroMagnitude = mag;
            }
            frequency[k] = k / (fftLength * spacingMm);
            ttf[k] = zeroMagnitude > 0 ? mag / zeroMagnitude : Double.NaN;
        }
        return new TtfResult(
            lsfRadiusMm,
            lsf,
            frequency,
            ttf,
            fftLength,
            frequency.length > 1 ? frequency[1] - frequency[0] : Double.NaN,
            1.0 / (2.0 * spacingMm),
            crossingFrequency(frequency, ttf, 0.5),
            crossingFrequency(frequency, ttf, 0.1)
        );
    }

    private double[] buildPointSpread64Style(double[] esf, double spacingMm, int length) {
        double[] psf = new double[length];
        if (esf == null || esf.length < 2 || length < 2) {
            return psf;
        }
        int last = Math.min(length - 2, esf.length - 1);
        boolean droppingEdge = esf[0] > esf[esf.length - 1];
        for (int i = 1; i <= last; i++) {
            double delta = esf[i] - esf[i - 1];
            psf[i] = (droppingEdge ? -delta : delta) / spacingMm;
        }
        psf[0] = 0.0;
        psf[length - 1] = 0.0;
        return psf;
    }

    private double crossingFrequency(double[] freq, double[] ttf, double target) {
        for (int i = 1; i < freq.length; i++) {
            double prev = ttf[i - 1];
            double cur = ttf[i];
            if (Double.isNaN(prev) || Double.isNaN(cur)) {
                continue;
            }
            if (prev >= target && cur <= target) {
                double denom = cur - prev;
                if (Math.abs(denom) < 1e-12) {
                    return freq[i];
                }
                double frac = (target - prev) / denom;
                return freq[i - 1] + frac * (freq[i] - freq[i - 1]);
            }
        }
        return Double.NaN;
    }

    private boolean exportAnalysisToCSV(int slice, double[] radiusPx, double[] radiusMm, double[] avgProfile, double[] segmentedProfile, TtfResult ttf) {
        SaveDialog sd = new SaveDialog("Save Circular Edge TTF CSV", String.format("CircularEdgeTTF_slice%03d", slice), ".csv");
        String fileName = sd.getFileName();
        String directory = sd.getDirectory();
        if (fileName == null || directory == null) {
            IJ.log("Circular edge TTF export cancelled.");
            clearStatusAndProgress();
            return false;
        }
        Path outPath = Paths.get(directory, fileName);
        StringBuilder sb = new StringBuilder();
        String lineSep = System.lineSeparator();
        sb.append("Radius_px,Radius_mm,AverageProfile_HU,SegmentedProfile_HU,LSF_Radius_mm,LSF_HU_per_mm,Frequency_cycles_per_mm,TTF")
          .append(lineSep);
        int rows = Math.max(Math.max(avgProfile.length, ttf.lsf.length), ttf.ttf.length);
        for (int i = 0; i < rows; i++) {
            appendValue(sb, i < radiusPx.length ? radiusPx[i] : Double.NaN);
            appendValue(sb, i < radiusMm.length ? radiusMm[i] : Double.NaN);
            appendValue(sb, i < avgProfile.length ? avgProfile[i] : Double.NaN);
            appendValue(sb, i < segmentedProfile.length ? segmentedProfile[i] : Double.NaN);
            appendValue(sb, i < ttf.lsfRadiusMm.length ? ttf.lsfRadiusMm[i] : Double.NaN);
            appendValue(sb, i < ttf.lsf.length ? ttf.lsf[i] : Double.NaN);
            appendValue(sb, i < ttf.frequency.length ? ttf.frequency[i] : Double.NaN);
            appendLastValue(sb, i < ttf.ttf.length ? ttf.ttf[i] : Double.NaN);
            sb.append(lineSep);
        }
        try {
            Files.write(outPath, sb.toString().getBytes(StandardCharsets.UTF_8));
            IJ.log("Circular edge TTF saved to " + outPath.toString());
            clearStatusAndProgress();
            return true;
        } catch (Exception e) {
            IJ.error("Export failed", e.getMessage());
            clearStatusAndProgress();
            return false;
        }
    }

    private boolean confirmFinishAfterCalculation(String analysisName, boolean csvSaved) {
        String csvStatus = csvSaved ? "CSVは保存済みです。" : "CSVは保存されていません。";
        YesNoCancelDialog dialog = new YesNoCancelDialog(
            IJ.getInstance(),
            analysisName + " calculation finished",
            "計算が完了しました。\n" + csvStatus + "\n" + analysisName + "をここで終了してよいですか？\n\n「続ける」またはキャンセルでは、次の計算を続けます。",
            "終了",
            "続ける"
        );
        return dialog.yesPressed();
    }

    protected void clearStatusAndProgress() {
        IJ.showProgress(1.0);
        IJ.showStatus("");
    }

    private void appendValue(StringBuilder sb, double value) {
        if (!Double.isNaN(value)) {
            sb.append(value);
        }
        sb.append(',');
    }

    private void appendLastValue(StringBuilder sb, double value) {
        if (!Double.isNaN(value)) {
            sb.append(value);
        }
    }

    private EdgeSegmentData extractSegmentData(SliceResult result) {
        int sIdx = (int) Math.ceil(result.highIntersection);
        int eIdx = (int) Math.floor(result.lowIntersection);
        sIdx = Math.max(0, Math.min(maxRadius, sIdx));
        eIdx = Math.max(0, Math.min(maxRadius, eIdx));
        List<Double> xbList = new ArrayList<>();
        List<Double> ybList = new ArrayList<>();
        xbList.add(0.0);
        ybList.add(result.highPlateau);
        xbList.add(result.highIntersection);
        ybList.add(result.highPlateau);
        if (result.lowIntersection > result.highIntersection && sIdx <= eIdx) {
            for (int r = sIdx; r <= eIdx; r++) {
                xbList.add((double) r);
                ybList.add(result.meanProfile[r]);
            }
        } else {
            int idx = (int) Math.round(Math.max(0, Math.min(maxRadius, result.highIntersection)));
            xbList.add(result.highIntersection);
            ybList.add(result.meanProfile[idx]);
        }
        xbList.add(result.lowIntersection);
        ybList.add(result.lowPlateau);
        xbList.add((double) maxRadius);
        ybList.add(result.lowPlateau);
        double[] xb = new double[xbList.size()];
        double[] yb = new double[ybList.size()];
        for (int i = 0; i < xb.length; i++) {
            xb[i] = xbList.get(i);
            yb[i] = ybList.get(i);
        }
        return new EdgeSegmentData(xb, yb);
    }

    private double interpolateProfile(double[] profile, double radiusPx) {
        if (profile == null || profile.length == 0) {
            return Double.NaN;
        }
        if (radiusPx <= 0) {
            return profile[0];
        }
        int last = profile.length - 1;
        if (radiusPx >= last) {
            return profile[last];
        }
        int r0 = (int) Math.floor(radiusPx);
        int r1 = r0 + 1;
        double frac = radiusPx - r0;
        return profile[r0] * (1.0 - frac) + profile[r1] * frac;
    }

    private float[] fitEdgeProfiles(FloatProcessor polar, double r0) {
        int w = polar.getWidth();
        int h = polar.getHeight();
        float[] offs = new float[w];
        double[] rawY = new double[h];
        double upperClipBase = (expectedEdgePx > 0 && tolerancePx > 0)
            ? Math.min(h - 1, expectedEdgePx + Math.max(tolerancePx, 3.0))
            : h - 1;
        double lowerClipBase = (expectedEdgePx > 0 && tolerancePx > 0)
            ? Math.max(0, expectedEdgePx - Math.max(tolerancePx, 3.0))
            : 0;

        for (int a = 0; a < w; a++) {
            List<Double> xList = new ArrayList<>();
            List<Double> yList = new ArrayList<>();
            for (int i = 0; i < h; i++) {
                double val = polar.getf(a, i);
                if (Double.isNaN(val) || val < OUTER_MIN_HU) {
                    val = OUTER_MIN_HU;
                }
                rawY[i] = val;
                if (i >= lowerClipBase && i <= upperClipBase) {
                    xList.add((double) i);
                    yList.add(val);
                }
            }

            double edgePos = Double.NaN;
            if (xList.size() >= 6) {
                double[] xFit = new double[xList.size()];
                double[] yFit = new double[yList.size()];
                for (int i = 0; i < xFit.length; i++) {
                    xFit[i] = xList.get(i);
                    yFit[i] = yList.get(i);
                }
                CurveFitter cf = new CurveFitter(xFit, yFit);
                double startLower = yFit[0];
                double startSlope = yFit[yFit.length - 1] - yFit[0];
                if (Math.abs(startSlope) < 1e-6) {
                    startSlope = 1.0;
                }
                cf.setInitialParameters(new double[]{startLower, startSlope, r0, 1.0});
                try {
                    cf.doFit(CurveFitter.RODBARD);
                    double[] params = cf.getParams();
                    edgePos = params[2];
                    if (Double.isNaN(edgePos) || edgePos < lowerClipBase || edgePos > upperClipBase) {
                        edgePos = Double.NaN;
                    }
                } catch (Exception e) {
                    edgePos = Double.NaN;
                }
            }

            if (Double.isNaN(edgePos)) {
                edgePos = estimateEdgeByGradient(rawY);
            }

            if (expectedEdgePx > 0 && tolerancePx > 0 && !Double.isNaN(edgePos)) {
                double clampLower = Math.max(0, expectedEdgePx - tolerancePx);
                double clampUpper = Math.min(h - 1, expectedEdgePx + tolerancePx);
                if (edgePos < clampLower) {
                    edgePos = clampLower;
                } else if (edgePos > clampUpper) {
                    edgePos = clampUpper;
                }
            }

            if (Double.isNaN(edgePos)) {
                edgePos = r0;
            }

            offs[a] = (float) (edgePos - r0);
        }
        return offs;
    }

    private float estimateEdgeByGradient(double[] profile) {
        int n = profile.length;
        if (n < 2) {
            return Float.NaN;
        }
        boolean dropping = profile[0] > profile[n - 1];
        int bestIdx = 0;
        double bestDiff = dropping ? Double.POSITIVE_INFINITY : Double.NEGATIVE_INFINITY;
        for (int i = 0; i < n - 1; i++) {
            double diff = profile[i + 1] - profile[i];
            if (dropping) {
                if (diff < bestDiff) {
                    bestDiff = diff;
                    bestIdx = i;
                }
            } else {
                if (diff > bestDiff) {
                    bestDiff = diff;
                    bestIdx = i;
                }
            }
        }
        double inner = profile[0];
        double outer = profile[n - 1];
        double target = (inner + outer) * 0.5;
        int start = Math.max(0, bestIdx - 3);
        int end = Math.min(n - 2, bestIdx + 3);
        for (int i = start; i <= end; i++) {
            double v1 = profile[i];
            double v2 = profile[i + 1];
            if (dropping ? (v1 >= target && v2 <= target) : (v1 <= target && v2 >= target)) {
                double denom = v2 - v1;
                if (Math.abs(denom) < 1e-6) {
                    return (float) i;
                }
                double frac = (target - v1) / denom;
                frac = Math.max(0.0, Math.min(1.0, frac));
                return (float) (i + frac);
            }
        }
        return (float) Math.max(0, Math.min(n - 1, bestIdx + 0.5));
    }

    private FloatProcessor alignPolarProfiles(FloatProcessor p, float[] offs) {
        int w = p.getWidth();
        int h = p.getHeight();
        FloatProcessor res = new FloatProcessor(w, h);
        for (int x = 0; x < w; x++) {
            for (int y = 0; y < h; y++) {
                double pos = y + offs[x];
                if (pos < 0) pos = 0;
                if (pos > h - 1) pos = h - 1;
                int y0 = (int) Math.floor(pos);
                int y1 = (int) Math.ceil(pos);
                float v0 = p.getf(x, y0);
                float v1 = p.getf(x, y1);
                double f = pos - y0;
                res.setf(x, y, (float) (v0 * (1 - f) + v1 * f));
            }
        }
        return res;
    }

    private static class EdgeSegmentData {
        final double[] radius;
        final double[] intensity;

        EdgeSegmentData(double[] radius, double[] intensity) {
            this.radius = radius;
            this.intensity = intensity;
        }
    }

    private static class PlotLimits {
        final double xMin;
        final double xMax;
        final double yMin;
        final double yMax;

        PlotLimits(double xMin, double xMax, double yMin, double yMax) {
            this.xMin = xMin;
            this.xMax = xMax;
            this.yMin = yMin;
            this.yMax = yMax;
        }
    }

    private class TtfCenterPreviewController implements TextListener {
        private final ImagePlus imp;
        private final TextField centerXField;
        private final TextField centerYField;
        private final TextField diameterMmField;
        private final TextField samplingRadiusField;
        private final Overlay originalOverlay;

        TtfCenterPreviewController(ImagePlus imp, TextField centerXField, TextField centerYField,
                                   TextField diameterMmField, TextField samplingRadiusField) {
            this.imp = imp;
            this.centerXField = centerXField;
            this.centerYField = centerYField;
            this.diameterMmField = diameterMmField;
            this.samplingRadiusField = samplingRadiusField;
            this.originalOverlay = imp.getOverlay();
        }

        void install() {
            centerXField.addTextListener(this);
            centerYField.addTextListener(this);
            diameterMmField.addTextListener(this);
            samplingRadiusField.addTextListener(this);
            updateFromFields();
        }

        void uninstall() {
            centerXField.removeTextListener(this);
            centerYField.removeTextListener(this);
            diameterMmField.removeTextListener(this);
            samplingRadiusField.removeTextListener(this);
            imp.setOverlay(originalOverlay);
            imp.updateAndDraw();
        }

        @Override
        public void textValueChanged(TextEvent e) {
            updateFromFields();
        }

        void updateFromFields() {
            PreviewGeometry geometry = currentGeometry();
            if (geometry == null) {
                imp.setOverlay(originalOverlay);
                imp.updateAndDraw();
                return;
            }

            boolean inside = geometry.centerX >= 0
                && geometry.centerX < imp.getWidth()
                && geometry.centerY >= 0
                && geometry.centerY < imp.getHeight()
                && geometry.samplingBounds.x >= 0
                && geometry.samplingBounds.y >= 0
                && geometry.samplingBounds.x + geometry.samplingBounds.width <= imp.getWidth()
                && geometry.samplingBounds.y + geometry.samplingBounds.height <= imp.getHeight();

            Overlay overlay = originalOverlay == null ? new Overlay() : originalOverlay.duplicate();
            addThumbnailPreview(overlay, geometry, inside);

            OvalRoi samplingRoi = new OvalRoi(
                geometry.centerX - geometry.samplingRadiusPx,
                geometry.centerY - geometry.samplingRadiusPx,
                geometry.samplingRadiusPx * 2.0,
                geometry.samplingRadiusPx * 2.0
            );
            samplingRoi.setName("TTF sampling radius preview");
            samplingRoi.setStrokeColor(inside ? Color.CYAN : Color.RED);
            samplingRoi.setStrokeWidth(2.0);
            overlay.add(samplingRoi);

            OvalRoi edgeRoi = new OvalRoi(
                geometry.centerX - geometry.edgeRadiusPx,
                geometry.centerY - geometry.edgeRadiusPx,
                geometry.edgeRadiusPx * 2.0,
                geometry.edgeRadiusPx * 2.0
            );
            edgeRoi.setName("TTF circular edge preview");
            edgeRoi.setStrokeColor(Color.YELLOW);
            edgeRoi.setStrokeWidth(2.0);
            overlay.add(edgeRoi);

            Roi centerMarker = new Roi(
                (int) Math.round(geometry.centerX) - 2,
                (int) Math.round(geometry.centerY) - 2,
                5,
                5
            );
            centerMarker.setName("TTF center preview");
            centerMarker.setStrokeColor(inside ? Color.GREEN : Color.RED);
            centerMarker.setStrokeWidth(2.0);
            overlay.add(centerMarker);

            imp.setOverlay(overlay);
            imp.updateAndDraw();
        }

        private PreviewGeometry currentGeometry() {
            try {
                double cx = Double.parseDouble(centerXField.getText());
                double cy = Double.parseDouble(centerYField.getText());
                double diameterMm = Double.parseDouble(diameterMmField.getText());
                double samplingRadius = Double.parseDouble(samplingRadiusField.getText());
                if (diameterMm <= 0 || samplingRadius <= 0 || pixelSizeMm <= 0) {
                    return null;
                }
                double edgeRadius = (diameterMm / 2.0) / pixelSizeMm;
                int x0 = (int) Math.round(cx - samplingRadius);
                int y0 = (int) Math.round(cy - samplingRadius);
                int size = Math.max(1, (int) Math.round(samplingRadius * 2.0));
                return new PreviewGeometry(cx, cy, edgeRadius, samplingRadius, new Rectangle(x0, y0, size, size));
            } catch (NumberFormatException e) {
                return null;
            }
        }

        private void addThumbnailPreview(Overlay overlay, PreviewGeometry geometry, boolean inside) {
            Rectangle imageBounds = new Rectangle(0, 0, imp.getWidth(), imp.getHeight());
            Rectangle cropBounds = geometry.samplingBounds.intersection(imageBounds);
            if (cropBounds.isEmpty()) {
                return;
            }

            int thumbnailSize = Math.max(32, Math.min(128, Math.min(imp.getWidth(), imp.getHeight()) / 5));
            Rectangle thumbnailBounds = thumbnailBounds(thumbnailSize, geometry.samplingBounds);
            try {
                ImageProcessor processor = imp.getProcessor().duplicate();
                processor.setMinAndMax(imp.getDisplayRangeMin(), imp.getDisplayRangeMax());
                processor.setRoi(cropBounds);
                ImageProcessor thumbnail = processor.crop().resize(thumbnailSize, thumbnailSize, true).convertToByte(true).convertToRGB();
                ImageRoi thumbnailRoi = new ImageRoi(thumbnailBounds.x, thumbnailBounds.y, thumbnail);
                thumbnailRoi.setName("TTF center thumbnail");
                overlay.add(thumbnailRoi);

                Roi border = new Roi(thumbnailBounds.x, thumbnailBounds.y, thumbnailBounds.width, thumbnailBounds.height);
                border.setName("TTF center thumbnail border");
                border.setStrokeColor(inside ? Color.CYAN : Color.RED);
                border.setStrokeWidth(2.0);
                overlay.add(border);
            } catch (RuntimeException e) {
                // Keep the circle preview even if thumbnail conversion fails in a display-specific image.
            }
        }

        private Rectangle thumbnailBounds(int thumbnailSize, Rectangle selectedArea) {
            int margin = Math.max(4, thumbnailSize / 12);
            Rectangle[] candidates = new Rectangle[] {
                new Rectangle(imp.getWidth() - thumbnailSize - margin, margin, thumbnailSize, thumbnailSize),
                new Rectangle(margin, margin, thumbnailSize, thumbnailSize),
                new Rectangle(imp.getWidth() - thumbnailSize - margin, imp.getHeight() - thumbnailSize - margin, thumbnailSize, thumbnailSize),
                new Rectangle(margin, imp.getHeight() - thumbnailSize - margin, thumbnailSize, thumbnailSize)
            };
            for (Rectangle candidate : candidates) {
                if (!candidate.intersects(selectedArea)) {
                    return candidate;
                }
            }
            return candidates[0];
        }
    }

    private static class PreviewGeometry {
        final double centerX;
        final double centerY;
        final double edgeRadiusPx;
        final double samplingRadiusPx;
        final Rectangle samplingBounds;

        PreviewGeometry(double centerX, double centerY, double edgeRadiusPx, double samplingRadiusPx, Rectangle samplingBounds) {
            this.centerX = centerX;
            this.centerY = centerY;
            this.edgeRadiusPx = edgeRadiusPx;
            this.samplingRadiusPx = samplingRadiusPx;
            this.samplingBounds = samplingBounds;
        }
    }

    private static class TtfImageCoordinateProbe extends MouseAdapter {
        private final ImagePlus imp;
        private final ImageCanvas canvas;
        private final Calibration cal;
        private final TextField centerXField;
        private final TextField centerYField;
        private final TtfCenterPreviewController preview;

        TtfImageCoordinateProbe(ImagePlus imp, ImageCanvas canvas, Calibration cal,
                                TextField centerXField, TextField centerYField,
                                TtfCenterPreviewController preview) {
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
            updateStatus(e, false);
        }

        @Override
        public void mouseClicked(MouseEvent e) {
            updateStatus(e, true);
        }

        private void updateStatus(MouseEvent e, boolean setCenter) {
            int x = canvas.offScreenX(e.getX());
            int y = canvas.offScreenY(e.getY());
            if (x < 0 || y < 0 || x >= imp.getWidth() || y >= imp.getHeight()) {
                return;
            }
            if (setCenter) {
                centerXField.setText(Integer.toString(x));
                centerYField.setText(Integer.toString(y));
                if (preview != null) {
                    preview.updateFromFields();
                }
            }
            double raw = imp.getProcessor().getf(x, y);
            double value = cal != null ? cal.getCValue(raw) : raw;
            IJ.showStatus(String.format(
                Locale.US,
                "Circular Edge center probe: x=%d, y=%d, value=%.3f%s",
                x,
                y,
                value,
                setCenter ? " -> center" : ""
            ));
        }
    }

    private class ProfileAdjustmentController {
        private final double[] axis;
        private final TextField highField;
        private final TextField lowField;
        private final TextField highIntersectionField;
        private final TextField lowIntersectionField;
        private SliceResult current;
        private PlotWindow[] windows;
        private ImageCanvas activeCanvas;
        private ProfileProbe activeProbe;

        ProfileAdjustmentController(double[] axis, SliceResult current, PlotWindow[] windows,
                                    TextField highField, TextField lowField,
                                    TextField highIntersectionField, TextField lowIntersectionField) {
            this.axis = axis;
            this.current = current;
            this.windows = windows;
            this.highField = highField;
            this.lowField = lowField;
            this.highIntersectionField = highIntersectionField;
            this.lowIntersectionField = lowIntersectionField;
        }

        void install() {
            installProbeOnFlattenedPlot();
        }

        PlotWindow[] getWindows() {
            return windows;
        }

        SliceResult getCurrent() {
            return current;
        }

        void redrawFromFields() {
            updateFromFieldsAndRedraw();
        }

        void applyClickedPlateau(ProfileSample sample, boolean useHigh) {
            TextField plateauTarget = useHigh ? highField : lowField;
            TextField intersectionTarget = useHigh ? highIntersectionField : lowIntersectionField;
            plateauTarget.setText(String.format(Locale.US, "%.3f", sample.valueHu));
            intersectionTarget.setText(String.format(Locale.US, "%.3f", sample.radiusPx));
            updateFromFieldsAndRedraw(false);
        }

        private void updateFromFieldsAndRedraw() {
            updateFromFieldsAndRedraw(true);
        }

        private void updateFromFieldsAndRedraw(boolean autoRecomputeIntersections) {
            double highPlateau = parseField(highField, current.highPlateau);
            double lowPlateau = parseField(lowField, current.lowPlateau);
            double highIntersection = parseField(highIntersectionField, current.highIntersection);
            double lowIntersection = parseField(lowIntersectionField, current.lowIntersection);

            boolean highPlateauChanged = Math.abs(highPlateau - current.highPlateau) > 1e-9;
            boolean lowPlateauChanged = Math.abs(lowPlateau - current.lowPlateau) > 1e-9;
            boolean highIntersectionChanged = Math.abs(highIntersection - current.highIntersection) > 1e-9;
            boolean lowIntersectionChanged = Math.abs(lowIntersection - current.lowIntersection) > 1e-9;

            if (autoRecomputeIntersections
                && ((highPlateauChanged && !highIntersectionChanged) || (lowPlateauChanged && !lowIntersectionChanged))) {
                IntersectionPair auto = recomputeIntersections(
                    current.meanProfile,
                    highPlateau,
                    lowPlateau,
                    current.highIntersection,
                    current.lowIntersection
                );
                if (highPlateauChanged && !highIntersectionChanged) {
                    highIntersection = auto.high;
                }
                if (lowPlateauChanged && !lowIntersectionChanged) {
                    lowIntersection = auto.low;
                }
            }

            highIntersection = Math.max(0, Math.min(maxRadius, highIntersection));
            lowIntersection = Math.max(0, Math.min(maxRadius, lowIntersection));
            if (lowIntersection < highIntersection) {
                double tmp = highIntersection;
                highIntersection = lowIntersection;
                lowIntersection = tmp;
            }

            highField.setText(String.format(Locale.US, "%.3f", highPlateau));
            lowField.setText(String.format(Locale.US, "%.3f", lowPlateau));
            highIntersectionField.setText(String.format(Locale.US, "%.3f", highIntersection));
            lowIntersectionField.setText(String.format(Locale.US, "%.3f", lowIntersection));

            current = new SliceResult(
                current.rawPolar,
                current.alignedPolar,
                current.meanProfile,
                highPlateau,
                lowPlateau,
                highIntersection,
                lowIntersection
            );
            removeActiveProbe();
            windows = showProfilePlots(current, axis, windows);
            installProbeOnFlattenedPlot();
        }

        private double parseField(TextField field, double fallback) {
            try {
                double value = Double.parseDouble(field.getText());
                return Double.isNaN(value) ? fallback : value;
            } catch (NumberFormatException e) {
                return fallback;
            }
        }

        private void installProbeOnFlattenedPlot() {
            if (windows == null || windows.length == 0 || windows[0] == null) {
                return;
            }
            ImageCanvas canvas = windows[0].getCanvas();
            Plot plot = windows[0].getPlot();
            if (canvas == null || plot == null) {
                return;
            }
            activeCanvas = canvas;
            activeProbe = new ProfileProbe(plot, canvas, current, this);
            canvas.addMouseMotionListener(activeProbe);
            canvas.addMouseListener(activeProbe);
        }

        private void removeActiveProbe() {
            if (activeCanvas != null && activeProbe != null) {
                activeCanvas.removeMouseMotionListener(activeProbe);
                activeCanvas.removeMouseListener(activeProbe);
            }
            activeCanvas = null;
            activeProbe = null;
        }
    }

    private class ProfileProbe extends MouseAdapter {
        private final Plot plot;
        private final ImageCanvas canvas;
        private final SliceResult result;
        private final ProfileAdjustmentController controller;
        private final double boundaryRadius;

        ProfileProbe(Plot plot, ImageCanvas canvas, SliceResult result, ProfileAdjustmentController controller) {
            this.plot = plot;
            this.canvas = canvas;
            this.result = result;
            this.controller = controller;
            this.boundaryRadius = expectedEdgePx > 0
                ? expectedEdgePx
                : (result.highIntersection + result.lowIntersection) / 2.0;
        }

        @Override
        public void mouseMoved(MouseEvent e) {
            ProfileSample sample = sample(e);
            if (sample != null) {
                IJ.showStatus(String.format(
                    Locale.US,
                    "Circular Edge probe: radius=%.2f px, curve=%.2f HU",
                    sample.radiusPx,
                    sample.valueHu
                ));
            }
        }

        @Override
        public void mouseClicked(MouseEvent e) {
            if (e.isPopupTrigger() || e.getButton() != MouseEvent.BUTTON1) {
                return;
            }
            ProfileSample sample = sample(e);
            if (sample == null || Double.isNaN(sample.valueHu)) {
                return;
            }

            boolean forceHigh = e.isShiftDown();
            boolean forceLow = e.isControlDown() || e.isAltDown();
            boolean useHigh = forceHigh || (!forceLow && sample.radiusPx <= boundaryRadius);
            String targetName = useHigh ? "High Plateau" : "Low Plateau";
            controller.applyClickedPlateau(sample, useHigh);

            String message = String.format(
                Locale.US,
                "Circular Edge probe: radius=%.2f px, curve=%.3f HU -> %s and intersection; plots redrawn",
                sample.radiusPx,
                sample.valueHu,
                targetName
            );
            IJ.showStatus(message);
            IJ.log(message);
        }

        private ProfileSample sample(MouseEvent e) {
            if (plot == null || canvas == null || result == null || result.meanProfile == null) {
                return null;
            }
            int mouseX = canvas.offScreenX(e.getX());
            double radius = plot.descaleX(mouseX);
            if (Double.isNaN(radius) || Double.isInfinite(radius)) {
                return null;
            }
            int last = Math.min(maxRadius, result.meanProfile.length - 1);
            radius = Math.max(0.0, Math.min(last, radius));
            double value = interpolateProfile(result.meanProfile, radius);
            return new ProfileSample(radius, value);
        }
    }

    private static class ProfileSample {
        final double radiusPx;
        final double valueHu;

        ProfileSample(double radiusPx, double valueHu) {
            this.radiusPx = radiusPx;
            this.valueHu = valueHu;
        }
    }

    private static class IntersectionPair {
        final double high;
        final double low;

        IntersectionPair(double high, double low) {
            this.high = high;
            this.low = low;
        }
    }

    protected static class PixelGeometry {
        final double pixelSizeMm;
        final double fovMm;
        final int matrix;
        final boolean isDicomDerived;
        final String source;

        PixelGeometry(double pixelSizeMm, double fovMm, int matrix, boolean isDicomDerived, String source) {
            this.pixelSizeMm = pixelSizeMm;
            this.fovMm = fovMm;
            this.matrix = matrix;
            this.isDicomDerived = isDicomDerived;
            this.source = source;
        }
    }

    protected static class SliceResult {
        final FloatProcessor rawPolar;
        final FloatProcessor alignedPolar;
        final double[] meanProfile;
        final double highPlateau;
        final double lowPlateau;
        final double highIntersection;
        final double lowIntersection;

        SliceResult(FloatProcessor rawPolar, FloatProcessor alignedPolar, double[] meanProfile,
                    double highPlateau, double lowPlateau, double highIntersection, double lowIntersection) {
            this.rawPolar = rawPolar;
            this.alignedPolar = alignedPolar;
            this.meanProfile = meanProfile;
            this.highPlateau = highPlateau;
            this.lowPlateau = lowPlateau;
            this.highIntersection = highIntersection;
            this.lowIntersection = lowIntersection;
        }
    }

    protected static class TtfResult {
        final double[] lsfRadiusMm;
        final double[] lsf;
        final double[] frequency;
        final double[] ttf;
        final int fftLength;
        final double frequencyStep;
        final double nyquist;
        final double ttf50;
        final double ttf10;

        TtfResult(double[] lsfRadiusMm, double[] lsf, double[] frequency, double[] ttf,
                  int fftLength, double frequencyStep, double nyquist, double ttf50, double ttf10) {
            this.lsfRadiusMm = lsfRadiusMm;
            this.lsf = lsf;
            this.frequency = frequency;
            this.ttf = ttf;
            this.fftLength = fftLength;
            this.frequencyStep = frequencyStep;
            this.nyquist = nyquist;
            this.ttf50 = ttf50;
            this.ttf10 = ttf10;
        }
    }
}
