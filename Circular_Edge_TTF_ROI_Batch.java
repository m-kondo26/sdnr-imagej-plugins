import ij.IJ;
import ij.ImagePlus;
import ij.gui.GenericDialog;
import ij.gui.NonBlockingGenericDialog;
import ij.gui.Overlay;
import ij.gui.OvalRoi;
import ij.gui.Plot;
import ij.gui.PlotWindow;
import ij.gui.Roi;
import ij.gui.YesNoCancelDialog;
import ij.io.OpenDialog;
import ij.io.RoiDecoder;
import ij.io.SaveDialog;
import ij.measure.Calibration;
import ij.process.ColorProcessor;
import ij.process.ImageProcessor;

import java.awt.Button;
import java.awt.Color;
import java.awt.event.ActionEvent;
import java.awt.geom.Rectangle2D;
import java.io.ByteArrayOutputStream;
import java.io.FileInputStream;
import java.io.IOException;
import java.lang.reflect.Field;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.List;
import java.util.Locale;
import java.util.zip.ZipEntry;
import java.util.zip.ZipInputStream;
import java.util.zip.ZipOutputStream;

/**
 * Batch circular-edge TTF/CTF measurement driven by named ImageJ ROI zip entries.
 */
public class Circular_Edge_TTF_ROI_Batch extends Circular_Edge_TTF {
    private static final String VERSION = "2026-09-09-negative-hu-fix";
    private static final int DEFAULT_FFT_SAMPLES = 64;
    private boolean autoRemainingButtonArmed = false;
    private boolean autoRemainingRequested = false;

    @Override
    public void run(String arg) {
        IJ.log("Circular Edge TTF ROI Batch " + VERSION + " loaded as class " + getClass().getName());
        ImagePlus imp = IJ.getImage();
        if (imp == null) {
            IJ.error("No image", "Open an image stack first.");
            return;
        }
        if (imp.getStackSize() < 1) {
            IJ.error("No stack", "The current image has no slices.");
            return;
        }

        Calibration cal = imp.getCalibration();
        PixelGeometry geometry = inferPixelGeometry(imp, cal, imp.getWidth(), imp.getHeight());
        if (!geometry.isDicomDerived) {
            IJ.log("Circular Edge TTF ROI Batch: DICOM FOV/pixel-spacing tags were not found. Source attempted: " + geometry.source);
            geometry = requestManualGeometry(geometry, imp.getWidth());
            if (geometry == null) {
                IJ.log("Circular Edge TTF ROI Batch cancelled: geometry input was cancelled.");
                return;
            }
        }
        applyGeometry(geometry);

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

        RoiEntry selected = chooseTargetRoi(roiEntries);
        if (selected == null) {
            return;
        }

        BatchParameters params = readBatchParameters(imp, selected, roiZipPath);
        if (params == null) {
            return;
        }
        applyBatchParameters(params);
        showSelectedRoiPreview(imp, params);

        List<SliceMeasurement> measurements = processSlices(imp, cal, params);
        if (measurements.isEmpty()) {
            IJ.log("Circular Edge TTF ROI Batch: no slices were measured.");
            clearStatusAndProgress();
            return;
        }

        BatchCurves curves = aggregateCurves(measurements);
        closeWindowByTitle("Polar Raw (HU)");
        closeWindowByTitle("Aligned Polar (HU)");
        closePlots(params.lastSliceWindows);
        showAggregatePlots(curves);
        boolean saved = params.exportWorkbook && exportBatchWorkbook(params, curves);
        IJ.log(String.format(
            Locale.US,
            "Circular Edge TTF ROI Batch finished: sample=%s slices=%d export=%s",
            params.roi.name,
            measurements.size(),
            saved ? "saved" : "not saved"
        ));
        clearStatusAndProgress();
    }

    private void applyGeometry(PixelGeometry geometry) {
        pixelSizeMm = geometry.pixelSizeMm;
        dicomFovMm = geometry.fovMm;
        dicomMatrix = geometry.matrix;
        geometryDicomDerived = geometry.isDicomDerived;
        pixelSizeSource = geometry.source;
        IJ.log(String.format(
            Locale.US,
            "Circular Edge TTF ROI Batch geometry: pixelSize=%.9f mm/pixel, FOV=%.6f mm, matrix=%d, dicomDerived=%s, source=%s",
            pixelSizeMm,
            dicomFovMm,
            dicomMatrix,
            geometryDicomDerived,
            pixelSizeSource
        ));
    }

    @Override
    protected void beforePlateauAdjustmentDialogShown(GenericDialog gd, SliceResult current, double[] axis) {
        if (!autoRemainingButtonArmed) {
            return;
        }
        gd.addButton("Use these positions for remaining slices", e -> {
            autoRemainingRequested = true;
            pressMeasureTtfButton(gd);
        });
    }

    @Override
    protected boolean acceptPlateauAdjustmentDialogDismissal() {
        return autoRemainingRequested;
    }

    private void pressMeasureTtfButton(GenericDialog gd) {
        Button ok = findMeasureTtfButton(gd);
        if (ok != null) {
            gd.actionPerformed(new ActionEvent(ok, ActionEvent.ACTION_PERFORMED, ok.getActionCommand()));
            return;
        }
        IJ.log("Circular Edge ROI Batch warning: Measure TTF button was not found; closing adjustment dialog directly.");
        gd.dispose();
    }

    private Button findMeasureTtfButton(GenericDialog gd) {
        Button[] buttons = gd.getButtons();
        if (buttons != null) {
            for (Button button : buttons) {
                if (button != null && isMeasureTtfButton(button)) {
                    return button;
                }
            }
        }
        try {
            Field okay = GenericDialog.class.getDeclaredField("okay");
            okay.setAccessible(true);
            Object value = okay.get(gd);
            if (value instanceof Button) {
                return (Button) value;
            }
        } catch (Exception e) {
            IJ.log("Circular Edge ROI Batch warning: could not access GenericDialog OK button: " + e.getMessage());
        }
        return null;
    }

    private boolean isMeasureTtfButton(Button button) {
        String label = button.getLabel();
        return "Measure TTF".equals(label) || "OK".equals(label);
    }

    private String chooseRoiZip() {
        OpenDialog od = new OpenDialog("Open ImageJ ROI zip", null);
        String fileName = od.getFileName();
        String directory = od.getDirectory();
        if (fileName == null || directory == null) {
            IJ.log("Circular Edge TTF ROI Batch cancelled: ROI zip was not selected.");
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
                    IJ.log("Circular Edge TTF ROI Batch: skipped unreadable ROI entry " + entry.getName());
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

    private RoiEntry chooseTargetRoi(List<RoiEntry> entries) {
        String[] names = new String[entries.size()];
        for (int i = 0; i < entries.size(); i++) {
            names[i] = entries.get(i).displayName;
        }
        GenericDialog gd = new GenericDialog("Select ROI Sample");
        gd.addMessage("Select the target sample ROI. The ROI centroid is used as the circular edge center.");
        gd.addChoice("Target sample ROI name", names, names[0]);
        gd.showDialog();
        if (gd.wasCanceled()) {
            return null;
        }
        String selected = gd.getNextChoice();
        for (RoiEntry entry : entries) {
            if (entry.displayName.equals(selected)) {
                return entry;
            }
        }
        return entries.get(0);
    }

    private BatchParameters readBatchParameters(ImagePlus imp, RoiEntry roi, String roiZipPath) {
        double defaultDiameterMm = defaultObjectDiameterMm(roi);
        double defaultToleranceMm = DEFAULT_TOLERANCE_MM;
        int defaultMaxRadius = defaultSamplingRadiusPx(defaultDiameterMm, defaultToleranceMm);
        int currentSlice = Math.max(1, imp.getCurrentSlice());

        GenericDialog gd = new NonBlockingGenericDialog("Circular Edge ROI Batch Parameters");
        gd.addMessage(String.format(
            Locale.US,
            "Sample ROI: %s%nCenter: x=%.3f, y=%.3f px%nGeometry: FOV %.6f mm / matrix %d = %.9f mm/pixel%nSource: %s",
            roi.name,
            roi.centerX,
            roi.centerY,
            dicomFovMm,
            dicomMatrix,
            pixelSizeMm,
            pixelSizeSource
        ));
        gd.addNumericField("Object diameter (mm)", defaultDiameterMm, 3);
        gd.addNumericField("Edge tolerance (mm)", defaultToleranceMm, 3);
        gd.addNumericField("Sampling radius (pixels)", defaultMaxRadius, 0);
        gd.addNumericField("Angular samples", DEFAULT_ANGLE_SAMPLES, 0);
        gd.addNumericField("TTF FFT samples", DEFAULT_FFT_SAMPLES, 0);
        gd.addNumericField("Start slice", currentSlice, 0);
        gd.addNumericField("End slice", imp.getStackSize(), 0);
        gd.addNumericField("Slice step", 1, 0);
        gd.addCheckbox("Confirm/correct upper and lower plateaus for each slice", true);
        gd.addCheckbox("Show polar raw/aligned images for each slice", false);
        gd.addCheckbox("Pause after each TTF/CTF plot", false);
        gd.addCheckbox("Export Excel workbook (.xlsx)", true);
        gd.showDialog();
        if (gd.wasCanceled()) {
            return null;
        }

        double enteredDiameter = gd.getNextNumber();
        double enteredTolerance = gd.getNextNumber();
        double enteredRadius = gd.getNextNumber();
        double enteredAngles = gd.getNextNumber();
        double enteredFftSamples = gd.getNextNumber();
        double enteredStart = gd.getNextNumber();
        double enteredEnd = gd.getNextNumber();
        double enteredStep = gd.getNextNumber();
        boolean confirmPlateaus = gd.getNextBoolean();
        boolean showPolar = gd.getNextBoolean();
        boolean pauseAfterSlice = gd.getNextBoolean();
        boolean export = gd.getNextBoolean();

        if (gd.invalidNumber()
            || enteredDiameter <= 0 || enteredTolerance <= 0 || enteredRadius < 8
            || enteredAngles < 16 || enteredFftSamples < 8
            || enteredStart < 1 || enteredEnd < 1 || enteredStep < 1
            || Double.isNaN(enteredDiameter) || Double.isNaN(enteredTolerance)
            || Double.isNaN(enteredRadius) || Double.isNaN(enteredAngles)
            || Double.isNaN(enteredFftSamples) || Double.isNaN(enteredStart)
            || Double.isNaN(enteredEnd) || Double.isNaN(enteredStep)) {
            IJ.error("Invalid input", "Diameter, tolerance, sampling radius, FFT samples, and slice range must be valid positive values.");
            return null;
        }
        double edgeRadiusPx = (enteredDiameter / 2.0) / pixelSizeMm;
        double toleranceRadiusPx = enteredTolerance / pixelSizeMm;
        if (edgeRadiusPx + toleranceRadiusPx >= enteredRadius - 1.0) {
            IJ.error(
                "Invalid sampling radius",
                "Sampling radius must be larger than object radius plus edge tolerance."
            );
            return null;
        }
        if (roi.centerX - enteredRadius < 0 || roi.centerY - enteredRadius < 0
            || roi.centerX + enteredRadius >= imp.getWidth()
            || roi.centerY + enteredRadius >= imp.getHeight()) {
            IJ.error(
                "Invalid sampling radius",
                "The selected ROI center and sampling radius extend outside the image."
            );
            return null;
        }

        int startSlice = clampSlice((int) Math.round(enteredStart), imp.getStackSize());
        int endSlice = clampSlice((int) Math.round(enteredEnd), imp.getStackSize());
        if (endSlice < startSlice) {
            IJ.error("Invalid slice range", "End slice must be greater than or equal to start slice.");
            return null;
        }
        int step = Math.max(1, (int) Math.round(enteredStep));

        return new BatchParameters(
            roi,
            roiZipPath,
            enteredDiameter,
            enteredTolerance,
            Math.max(8, (int) Math.round(enteredRadius)),
            Math.max(16, (int) Math.round(enteredAngles)),
            Math.max(8, (int) Math.round(enteredFftSamples)),
            startSlice,
            endSlice,
            step,
            confirmPlateaus,
            showPolar,
            pauseAfterSlice,
            export
        );
    }

    private void applyBatchParameters(BatchParameters params) {
        centerX = params.roi.centerX;
        centerY = params.roi.centerY;
        objectDiameterMm = params.objectDiameterMm;
        toleranceMm = params.toleranceMm;
        maxRadius = params.maxRadius;
        angleSamples = params.angleSamples;
        fftSamples = params.fftSamples;
        expectedEdgePx = (objectDiameterMm / 2.0) / pixelSizeMm;
        tolerancePx = Math.max(1.0, toleranceMm / pixelSizeMm);
        IJ.log(String.format(
            Locale.US,
            "Circular Edge TTF ROI Batch parameters: sample=%s center=(%.3f,%.3f) diameter=%.6f mm tolerance=%.6f mm maxRadius=%d angles=%d fft=%d slices=%d-%d step=%d confirm=%s",
            params.roi.name,
            centerX,
            centerY,
            objectDiameterMm,
            toleranceMm,
            maxRadius,
            angleSamples,
            fftSamples,
            params.startSlice,
            params.endSlice,
            params.sliceStep,
            params.confirmPlateaus
        ));
    }

    private List<SliceMeasurement> processSlices(ImagePlus imp, Calibration cal, BatchParameters params) {
        List<SliceMeasurement> measurements = new ArrayList<>();
        int radialLength = maxRadius + 1;
        double[] radiusPx = buildRadiusAxis(radialLength);
        int planned = ((params.endSlice - params.startSlice) / params.sliceStep) + 1;
        int completed = 0;

        for (int slice = params.startSlice; slice <= params.endSlice; slice += params.sliceStep) {
            IJ.showProgress(completed, planned);
            IJ.showStatus("Circular Edge ROI Batch: measuring slice " + slice);
            imp.setSlice(slice);

            ImageProcessor ip = imp.getStack().getProcessor(slice).duplicate();
            if (imp.getType() == ImagePlus.COLOR_RGB || ip instanceof ColorProcessor) {
                ip = ip.convertToByte(true);
            }
            ip.setInterpolationMethod(ImageProcessor.BILINEAR);

            SliceResult autoResult;
            try {
                autoResult = analyzeSlice(ip, cal, imp.getWidth(), imp.getHeight(), radialLength);
            } catch (IllegalArgumentException e) {
                IJ.error("Edge alignment failed", "Slice " + slice + ": " + e.getMessage()
                    + "\nThe batch was stopped; no workbook will be exported from this run.");
                // Do not present an incomplete batch as a completed measurement.
                measurements.clear();
                break;
            }
            if (params.showPolarImages) {
                autoResult.rawPolar.resetMinAndMax();
                autoResult.alignedPolar.resetMinAndMax();
                showOrReplaceImage("Polar Raw (HU)", autoResult.rawPolar);
                showOrReplaceImage("Aligned Polar (HU)", autoResult.alignedPolar);
            }

            SliceResult tuned;
            if (params.autoRemaining) {
                tuned = resultFromFixedIntersections(
                    autoResult,
                    params.fixedHighIntersection,
                    params.fixedLowIntersection
                );
            } else if (params.confirmPlateaus) {
                SliceResult proposedResult = !measurements.isEmpty()
                    ? proposalFromPreviousSlice(autoResult, measurements.get(measurements.size() - 1))
                    : autoResult;
                autoRemainingRequested = false;
                autoRemainingButtonArmed = true;
                tuned = adjustPlateausInteractively(proposedResult, radiusPx);
                autoRemainingButtonArmed = false;
                if (tuned != null && autoRemainingRequested) {
                    params.autoRemaining = true;
                    params.fixedHighIntersection = tuned.highIntersection;
                    params.fixedLowIntersection = tuned.lowIntersection;
                    IJ.log(String.format(
                        Locale.US,
                        "Circular Edge ROI Batch: remaining slices will use fixed intersections HighInt=%.3f px, LowInt=%.3f px.",
                        params.fixedHighIntersection,
                        params.fixedLowIntersection
                    ));
                }
            } else {
                tuned = autoResult;
            }
            if (tuned == null) {
                IJ.log("Circular Edge TTF ROI Batch cancelled at slice " + slice + ".");
                break;
            }

            double[] segmented = buildSegmentedProfile(
                tuned.meanProfile,
                tuned.highPlateau,
                tuned.lowPlateau,
                tuned.highIntersection,
                tuned.lowIntersection
            );
            TtfResult ttf = computeTtf(segmented, pixelSizeMm);
            double contrastHu = Math.abs(tuned.highPlateau - tuned.lowPlateau);
            double[] ctf = computeCtf(ttf.ttf, contrastHu);
            SliceMeasurement measurement = new SliceMeasurement(slice, tuned, ttf, ctf, contrastHu);
            measurements.add(measurement);
            params.lastSliceWindows = showSliceTransferPlots(measurement, params.lastSliceWindows);

            IJ.log(String.format(
                Locale.US,
                "Circular Edge ROI Batch slice %d: High=%.3f HU Low=%.3f HU Contrast=%.3f HU HighInt=%.3f px LowInt=%.3f px",
                slice,
                tuned.highPlateau,
                tuned.lowPlateau,
                contrastHu,
                tuned.highIntersection,
                tuned.lowIntersection
            ));

            completed++;
            if (params.pauseAfterSlice && slice + params.sliceStep <= params.endSlice && !confirmContinue(slice)) {
                break;
            }
        }

        IJ.showProgress(1.0);
        return measurements;
    }

    private PlotWindow[] showSliceTransferPlots(SliceMeasurement measurement, PlotWindow[] previous) {
        closePlots(previous);

        TtfResult ttf = measurement.ttf;
        Plot ttfPlot = new Plot(
            "TTF slice " + measurement.slice,
            "Spatial frequency (cycles/mm)",
            "TTF"
        );
        ttfPlot.setFrameSize(720, 360);
        ttfPlot.setLimits(0, ttf.nyquist, 0, 1.1);
        ttfPlot.setColor(Color.BLUE);
        ttfPlot.add("line", ttf.frequency, ttf.ttf);
        ttfPlot.setColor(Color.GRAY);
        ttfPlot.add("line", new double[]{0, ttf.nyquist}, new double[]{0.5, 0.5});
        PlotWindow ttfWindow = ttfPlot.show();
        if (ttfWindow != null) {
            ttfWindow.setLocation(100, 120);
            ttfWindow.toFront();
        }

        Plot ctfPlot = new Plot(
            "CTF slice " + measurement.slice,
            "Spatial frequency (cycles/mm)",
            "CTF (HU)"
        );
        ctfPlot.setFrameSize(720, 360);
        double maxCtf = maxFinite(measurement.ctf);
        ctfPlot.setLimits(0, ttf.nyquist, 0, Math.max(1.0, maxCtf * 1.1));
        ctfPlot.setColor(Color.RED);
        ctfPlot.add("line", ttf.frequency, measurement.ctf);
        PlotWindow ctfWindow = ctfPlot.show();
        if (ctfWindow != null) {
            ctfWindow.setLocation(ttfWindow != null ? ttfWindow.getX() + ttfWindow.getWidth() + 40 : 500, 120);
            ctfWindow.toFront();
        }

        return new PlotWindow[]{ttfWindow, ctfWindow};
    }

    private SliceResult proposalFromPreviousSlice(SliceResult current, SliceMeasurement previous) {
        double highIntersection = clampRadius(previous.result.highIntersection);
        double lowIntersection = clampRadius(previous.result.lowIntersection);
        SliceResult proposed = resultFromFixedIntersections(current, highIntersection, lowIntersection);
        IJ.log(String.format(
            Locale.US,
            "Circular Edge ROI Batch slice proposal: using previous slice %d intersections as initial values: HighInt=%.3f px, LowInt=%.3f px",
            previous.slice,
            proposed.highIntersection,
            proposed.lowIntersection
        ));
        return proposed;
    }

    private SliceResult resultFromFixedIntersections(SliceResult current, double highIntersection, double lowIntersection) {
        highIntersection = clampRadius(highIntersection);
        lowIntersection = clampRadius(lowIntersection);
        if (lowIntersection < highIntersection) {
            double tmp = highIntersection;
            highIntersection = lowIntersection;
            lowIntersection = tmp;
        }

        double highPlateau = interpolateCurrentProfile(current.meanProfile, highIntersection, current.highPlateau);
        double lowPlateau = interpolateCurrentProfile(current.meanProfile, lowIntersection, current.lowPlateau);
        return new SliceResult(
            current.rawPolar,
            current.alignedPolar,
            current.meanProfile,
            highPlateau,
            lowPlateau,
            highIntersection,
            lowIntersection
        );
    }

    private double clampRadius(double radius) {
        if (Double.isNaN(radius) || Double.isInfinite(radius)) {
            return 0.0;
        }
        return Math.max(0.0, Math.min(maxRadius, radius));
    }

    private double interpolateCurrentProfile(double[] profile, double radiusPx, double fallback) {
        if (profile == null || profile.length == 0 || Double.isNaN(radiusPx) || Double.isInfinite(radiusPx)) {
            return fallback;
        }
        if (radiusPx <= 0) {
            return finiteOrFallback(profile[0], fallback);
        }
        int last = Math.min(maxRadius, profile.length - 1);
        if (radiusPx >= last) {
            return finiteOrFallback(profile[last], fallback);
        }
        int r0 = (int) Math.floor(radiusPx);
        int r1 = Math.min(last, r0 + 1);
        double frac = radiusPx - r0;
        double value = profile[r0] * (1.0 - frac) + profile[r1] * frac;
        return finiteOrFallback(value, fallback);
    }

    private double finiteOrFallback(double value, double fallback) {
        return Double.isFinite(value) ? value : fallback;
    }

    private BatchCurves aggregateCurves(List<SliceMeasurement> measurements) {
        double[] frequency = measurements.get(0).ttf.frequency;
        int bins = frequency.length;
        int slices = measurements.size();
        double[][] ttfPerSlice = new double[slices][bins];
        double[][] ctfPerSlice = new double[slices][bins];
        int[] sliceNumbers = new int[slices];

        for (int s = 0; s < slices; s++) {
            SliceMeasurement measurement = measurements.get(s);
            sliceNumbers[s] = measurement.slice;
            for (int i = 0; i < bins; i++) {
                ttfPerSlice[s][i] = i < measurement.ttf.ttf.length ? measurement.ttf.ttf[i] : Double.NaN;
                ctfPerSlice[s][i] = i < measurement.ctf.length ? measurement.ctf[i] : Double.NaN;
            }
        }

        double[] averageTtf = new double[bins];
        double[] averageCtf = new double[bins];
        double[] sdTtf = new double[bins];
        double[] sdCtf = new double[bins];
        int[] nTtf = new int[bins];
        int[] nCtf = new int[bins];
        for (int i = 0; i < bins; i++) {
            averageTtf[i] = meanAt(ttfPerSlice, i);
            averageCtf[i] = meanAt(ctfPerSlice, i);
            sdTtf[i] = sdAt(ttfPerSlice, i, averageTtf[i]);
            sdCtf[i] = sdAt(ctfPerSlice, i, averageCtf[i]);
            nTtf[i] = countFiniteAt(ttfPerSlice, i);
            nCtf[i] = countFiniteAt(ctfPerSlice, i);
        }

        return new BatchCurves(
            frequency,
            ttfPerSlice,
            ctfPerSlice,
            averageTtf,
            averageCtf,
            sdTtf,
            sdCtf,
            nTtf,
            nCtf,
            sliceNumbers,
            measurements.get(0).ttf.nyquist
        );
    }

    private void showAggregatePlots(BatchCurves curves) {
        closeWindowByTitle("Circular Edge TTF by Slice");
        closeWindowByTitle("Circular Edge CTF by Slice");

        Plot ttfPlot = new Plot("Circular Edge TTF by Slice", "Spatial frequency (cycles/mm)", "TTF");
        ttfPlot.setFrameSize(820, 440);
        ttfPlot.setLimits(0, curves.nyquist, 0, 1.1);
        ttfPlot.setLineWidth(1);
        for (int s = 0; s < curves.ttfPerSlice.length; s++) {
            ttfPlot.setColor(sliceColor(s, curves.ttfPerSlice.length));
            ttfPlot.add("line", curves.frequency, curves.ttfPerSlice[s]);
        }
        ttfPlot.setLineWidth(3);
        ttfPlot.setColor(Color.BLACK);
        ttfPlot.add("line", curves.frequency, curves.averageTtf);
        ttfPlot.setLineWidth(1);
        PlotWindow ttfWindow = ttfPlot.show();
        if (ttfWindow != null) {
            ttfWindow.setLocation(80, 120);
            ttfWindow.toFront();
        }

        Plot ctfPlot = new Plot("Circular Edge CTF by Slice", "Spatial frequency (cycles/mm)", "CTF (HU)");
        ctfPlot.setFrameSize(820, 440);
        double maxCtf = maxFinite(curves.ctfPerSlice, curves.averageCtf);
        ctfPlot.setLimits(0, curves.nyquist, 0, Math.max(1.0, maxCtf * 1.1));
        ctfPlot.setLineWidth(1);
        for (int s = 0; s < curves.ctfPerSlice.length; s++) {
            ctfPlot.setColor(sliceColor(s, curves.ctfPerSlice.length));
            ctfPlot.add("line", curves.frequency, curves.ctfPerSlice[s]);
        }
        ctfPlot.setLineWidth(3);
        ctfPlot.setColor(Color.BLACK);
        ctfPlot.add("line", curves.frequency, curves.averageCtf);
        ctfPlot.setLineWidth(1);
        PlotWindow ctfWindow = ctfPlot.show();
        if (ctfWindow != null) {
            ctfWindow.setLocation(ttfWindow != null ? ttfWindow.getX() + ttfWindow.getWidth() + 40 : 940, 120);
            ctfWindow.toFront();
        }

        IJ.log("Circular Edge ROI Batch aggregate plots shown. Thin colored lines are per-slice curves; thick black line is the average.");
    }

    private boolean exportBatchWorkbook(BatchParameters params, BatchCurves curves) {
        SaveDialog sd = new SaveDialog(
            "Save Circular Edge TTF/CTF Excel workbook",
            "CircularEdgeTTF_CTF_" + safeFileStem(params.roi.name) + "_slices" + params.startSlice + "-" + params.endSlice,
            ".xlsx"
        );
        String fileName = sd.getFileName();
        String directory = sd.getDirectory();
        if (fileName == null || directory == null) {
            IJ.log("Circular Edge ROI Batch Excel export cancelled.");
            return false;
        }

        Path workbookPath = Paths.get(directory, ensureXlsxExtension(fileName));
        try {
            writeXlsxWorkbook(workbookPath, curves);
            IJ.log("Circular Edge ROI Batch Excel workbook saved to " + workbookPath.toString());
            return true;
        } catch (IOException e) {
            IJ.error("Export failed", e.getMessage());
            return false;
        }
    }

    private void writeXlsxWorkbook(Path path, BatchCurves curves) throws IOException {
        ZipOutputStream zip = new ZipOutputStream(Files.newOutputStream(path));
        try {
            putZipEntry(zip, "[Content_Types].xml", contentTypesXml());
            putZipEntry(zip, "_rels/.rels", rootRelationshipsXml());
            putZipEntry(zip, "xl/workbook.xml", workbookXml());
            putZipEntry(zip, "xl/_rels/workbook.xml.rels", workbookRelationshipsXml());
            putZipEntry(zip, "xl/worksheets/sheet1.xml", waveformSheetXml(
                "Average_TTF",
                curves.frequency,
                curves.averageTtf,
                curves.ttfPerSlice,
                curves.sliceNumbers
            ));
            putZipEntry(zip, "xl/worksheets/sheet2.xml", waveformSheetXml(
                "Average_CTF_HU",
                curves.frequency,
                curves.averageCtf,
                curves.ctfPerSlice,
                curves.sliceNumbers
            ));
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
            + "<Override PartName=\"/xl/worksheets/sheet2.xml\" ContentType=\"application/vnd.openxmlformats-officedocument.spreadsheetml.worksheet+xml\"/>"
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
            + "<sheet name=\"TTF\" sheetId=\"1\" r:id=\"rId1\"/>"
            + "<sheet name=\"CTF\" sheetId=\"2\" r:id=\"rId2\"/>"
            + "</sheets>"
            + "</workbook>";
    }

    private String workbookRelationshipsXml() {
        return "<?xml version=\"1.0\" encoding=\"UTF-8\"?>"
            + "<Relationships xmlns=\"http://schemas.openxmlformats.org/package/2006/relationships\">"
            + "<Relationship Id=\"rId1\" Type=\"http://schemas.openxmlformats.org/officeDocument/2006/relationships/worksheet\" Target=\"worksheets/sheet1.xml\"/>"
            + "<Relationship Id=\"rId2\" Type=\"http://schemas.openxmlformats.org/officeDocument/2006/relationships/worksheet\" Target=\"worksheets/sheet2.xml\"/>"
            + "</Relationships>";
    }

    private String waveformSheetXml(String averageHeader, double[] frequency, double[] average,
                                    double[][] perSlice, int[] sliceNumbers) {
        StringBuilder sb = new StringBuilder();
        sb.append("<?xml version=\"1.0\" encoding=\"UTF-8\"?>")
          .append("<worksheet xmlns=\"http://schemas.openxmlformats.org/spreadsheetml/2006/main\">")
          .append("<sheetData>");

        appendWaveformHeaderRow(sb, averageHeader, sliceNumbers);
        for (int i = 0; i < frequency.length; i++) {
            int row = i + 2;
            sb.append("<row r=\"").append(row).append("\">");
            appendNumberCell(sb, row, 1, frequency[i]);
            appendNumberCell(sb, row, 2, average[i]);
            for (int s = 0; s < sliceNumbers.length; s++) {
                double value = i < perSlice[s].length ? perSlice[s][i] : Double.NaN;
                appendNumberCell(sb, row, s + 3, value);
            }
            sb.append("</row>");
        }

        sb.append("</sheetData></worksheet>");
        return sb.toString();
    }

    private void appendWaveformHeaderRow(StringBuilder sb, String averageHeader, int[] sliceNumbers) {
        sb.append("<row r=\"1\">");
        appendInlineStringCell(sb, 1, 1, "Spatial_frequency_cycles_per_mm");
        appendInlineStringCell(sb, 1, 2, averageHeader);
        for (int s = 0; s < sliceNumbers.length; s++) {
            appendInlineStringCell(sb, 1, s + 3, "Slice_" + sliceNumbers[s]);
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
            } else {
                sb.append(ch);
            }
        }
        return sb.toString();
    }

    private String ensureXlsxExtension(String fileName) {
        return fileName.toLowerCase(Locale.US).endsWith(".xlsx") ? fileName : fileName + ".xlsx";
    }

    private boolean exportBatchCsv(BatchParameters params, List<SliceMeasurement> measurements, BatchCurves curves) {
        SaveDialog sd = new SaveDialog(
            "Save Circular Edge TTF/CTF CSV",
            "CircularEdgeTTF_CTF_" + safeFileStem(params.roi.name) + "_slices" + params.startSlice + "-" + params.endSlice,
            ".csv"
        );
        String fileName = sd.getFileName();
        String directory = sd.getDirectory();
        if (fileName == null || directory == null) {
            IJ.log("Circular Edge ROI Batch CSV export cancelled.");
            return false;
        }

        Path frequencyCsv = Paths.get(directory, fileName);
        Path summaryCsv = companionSummaryPath(frequencyCsv);
        try {
            writeUtf8BomCsv(frequencyCsv, buildFrequencyCsv(params, curves));
            writeUtf8BomCsv(summaryCsv, buildSummaryCsv(params, measurements));
            IJ.log("Circular Edge ROI Batch frequency CSV saved to " + frequencyCsv.toString());
            IJ.log("Circular Edge ROI Batch slice summary CSV saved to " + summaryCsv.toString());
            return true;
        } catch (IOException e) {
            IJ.error("Export failed", e.getMessage());
            return false;
        }
    }

    private String buildFrequencyCsv(BatchParameters params, BatchCurves curves) {
        StringBuilder sb = new StringBuilder();
        List<String> header = new ArrayList<>();
        header.add("Sample");
        header.add("Frequency_cycles_per_mm");
        header.add("TTF_Average");
        header.add("TTF_SD");
        header.add("TTF_N");
        header.add("CTF_Average_HU");
        header.add("CTF_SD_HU");
        header.add("CTF_N");
        for (int s = 0; s < curves.sliceNumbers.length; s++) {
            header.add("TTF_Slice_" + curves.sliceNumbers[s]);
            header.add("CTF_Slice_" + curves.sliceNumbers[s] + "_HU");
        }
        appendCsvRow(sb, header);

        for (int i = 0; i < curves.frequency.length; i++) {
            List<String> row = new ArrayList<>();
            row.add(params.roi.name);
            row.add(numberCell(curves.frequency[i]));
            row.add(numberCell(curves.averageTtf[i]));
            row.add(numberCell(curves.sdTtf[i]));
            row.add(Integer.toString(curves.nTtf[i]));
            row.add(numberCell(curves.averageCtf[i]));
            row.add(numberCell(curves.sdCtf[i]));
            row.add(Integer.toString(curves.nCtf[i]));
            for (int s = 0; s < curves.sliceNumbers.length; s++) {
                row.add(numberCell(curves.ttfPerSlice[s][i]));
                row.add(numberCell(curves.ctfPerSlice[s][i]));
            }
            appendCsvRow(sb, row);
        }
        return sb.toString();
    }

    private String buildSummaryCsv(BatchParameters params, List<SliceMeasurement> measurements) {
        StringBuilder sb = new StringBuilder();
        List<String> header = new ArrayList<>();
        header.add("Sample");
        header.add("ROI_Zip");
        header.add("ROI_Entry");
        header.add("Slice");
        header.add("Center_X_px");
        header.add("Center_Y_px");
        header.add("Object_Diameter_mm");
        header.add("Tolerance_mm");
        header.add("Sampling_Radius_px");
        header.add("Angular_Samples");
        header.add("TTF_FFT_N");
        header.add("PixelSize_mm");
        header.add("FOV_mm");
        header.add("Matrix");
        header.add("Geometry_DICOM_Derived");
        header.add("HighPlateau_HU");
        header.add("LowPlateau_HU");
        header.add("Contrast_HU_abs");
        header.add("HighIntersection_px");
        header.add("LowIntersection_px");
        header.add("TTF50_cycles_per_mm");
        header.add("TTF10_cycles_per_mm");
        appendCsvRow(sb, header);

        for (SliceMeasurement measurement : measurements) {
            List<String> row = new ArrayList<>();
            row.add(params.roi.name);
            row.add(params.roiZipPath);
            row.add(params.roi.entryName);
            row.add(Integer.toString(measurement.slice));
            row.add(numberCell(centerX));
            row.add(numberCell(centerY));
            row.add(numberCell(objectDiameterMm));
            row.add(numberCell(toleranceMm));
            row.add(Integer.toString(maxRadius));
            row.add(Integer.toString(angleSamples));
            row.add(Integer.toString(measurement.ttf.fftLength));
            row.add(numberCell(pixelSizeMm));
            row.add(numberCell(dicomFovMm));
            row.add(Integer.toString(dicomMatrix));
            row.add(geometryDicomDerived ? "1" : "0");
            row.add(numberCell(measurement.result.highPlateau));
            row.add(numberCell(measurement.result.lowPlateau));
            row.add(numberCell(measurement.contrastHu));
            row.add(numberCell(measurement.result.highIntersection));
            row.add(numberCell(measurement.result.lowIntersection));
            row.add(numberCell(measurement.ttf.ttf50));
            row.add(numberCell(measurement.ttf.ttf10));
            appendCsvRow(sb, row);
        }
        return sb.toString();
    }

    private void writeUtf8BomCsv(Path path, String text) throws IOException {
        Files.write(path, ("\uFEFF" + text).getBytes(StandardCharsets.UTF_8));
    }

    private double defaultObjectDiameterMm(RoiEntry roi) {
        double diameterPx = roi.diameterPx;
        if (diameterPx > 0 && pixelSizeMm > 0) {
            return diameterPx * pixelSizeMm;
        }
        return DEFAULT_DIAMETER_MM;
    }

    private int defaultSamplingRadiusPx(double diameterMm, double toleranceMm) {
        double edgeRadiusPx = (diameterMm / 2.0) / pixelSizeMm;
        double toleranceRadiusPx = toleranceMm / pixelSizeMm;
        return Math.max(DEFAULT_MAX_RADIUS, (int) Math.ceil(edgeRadiusPx + toleranceRadiusPx + 8.0));
    }

    private double[] buildRadiusAxis(int length) {
        double[] radius = new double[length];
        for (int i = 0; i < length; i++) {
            radius[i] = i;
        }
        return radius;
    }

    private double[] computeCtf(double[] ttf, double contrastHu) {
        double[] ctf = new double[ttf.length];
        for (int i = 0; i < ttf.length; i++) {
            ctf[i] = Double.isNaN(ttf[i]) ? Double.NaN : ttf[i] * contrastHu;
        }
        return ctf;
    }

    private boolean confirmContinue(int slice) {
        YesNoCancelDialog dialog = new YesNoCancelDialog(
            IJ.getInstance(),
            "Circular Edge ROI Batch",
            "Slice " + slice + " TTF/CTF plot is complete.\nContinue to the next slice?",
            "Continue",
            "Stop"
        );
        return dialog.yesPressed();
    }

    private void showSelectedRoiPreview(ImagePlus imp, BatchParameters params) {
        Overlay overlay = new Overlay();
        try {
            Roi selected = (Roi) params.roi.roi.clone();
            selected.setName("Selected ROI: " + params.roi.name);
            selected.setStrokeColor(Color.GREEN);
            selected.setStrokeWidth(2.0);
            overlay.add(selected);
        } catch (RuntimeException e) {
            IJ.log("Circular Edge ROI Batch: selected ROI preview skipped: " + e.getMessage());
        }

        double edgeRadius = (params.objectDiameterMm / 2.0) / pixelSizeMm;
        OvalRoi edge = new OvalRoi(
            params.roi.centerX - edgeRadius,
            params.roi.centerY - edgeRadius,
            edgeRadius * 2.0,
            edgeRadius * 2.0
        );
        edge.setName("Expected circular edge");
        edge.setStrokeColor(Color.YELLOW);
        edge.setStrokeWidth(2.0);
        overlay.add(edge);

        OvalRoi sampling = new OvalRoi(
            params.roi.centerX - params.maxRadius,
            params.roi.centerY - params.maxRadius,
            params.maxRadius * 2.0,
            params.maxRadius * 2.0
        );
        sampling.setName("TTF sampling radius");
        sampling.setStrokeColor(Color.CYAN);
        sampling.setStrokeWidth(1.0);
        overlay.add(sampling);
        imp.setOverlay(overlay);
        imp.updateAndDraw();
    }

    private int clampSlice(int slice, int stackSize) {
        return Math.max(1, Math.min(stackSize, slice));
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

    private String safeFileStem(String name) {
        String stem = name == null ? "ROI" : name.trim();
        if (stem.isEmpty()) {
            stem = "ROI";
        }
        stem = stem.replaceAll("[\\\\/:*?\"<>|]+", "_");
        stem = stem.replaceAll("\\s+", "_");
        return stem;
    }

    private Path companionSummaryPath(Path frequencyCsv) {
        String fileName = frequencyCsv.getFileName().toString();
        String lower = fileName.toLowerCase(Locale.US);
        String summaryName = lower.endsWith(".csv")
            ? fileName.substring(0, fileName.length() - 4) + "_slice_summary.csv"
            : fileName + "_slice_summary.csv";
        Path parent = frequencyCsv.getParent();
        return parent == null ? Paths.get(summaryName) : parent.resolve(summaryName);
    }

    private void appendCsvRow(StringBuilder sb, List<String> cells) {
        for (int i = 0; i < cells.size(); i++) {
            if (i > 0) {
                sb.append(',');
            }
            appendCsvCell(sb, cells.get(i));
        }
        sb.append(System.lineSeparator());
    }

    private void appendCsvCell(StringBuilder sb, String value) {
        if (value == null) {
            return;
        }
        boolean quote = value.indexOf(',') >= 0 || value.indexOf('"') >= 0
            || value.indexOf('\n') >= 0 || value.indexOf('\r') >= 0;
        if (!quote) {
            sb.append(value);
            return;
        }
        sb.append('"');
        for (int i = 0; i < value.length(); i++) {
            char ch = value.charAt(i);
            if (ch == '"') {
                sb.append("\"\"");
            } else {
                sb.append(ch);
            }
        }
        sb.append('"');
    }

    private String numberCell(double value) {
        return Double.isNaN(value) || Double.isInfinite(value)
            ? ""
            : String.format(Locale.US, "%.12g", value);
    }

    private double meanAt(double[][] values, int index) {
        double sum = 0.0;
        int n = 0;
        for (double[] row : values) {
            if (index < row.length && Double.isFinite(row[index])) {
                sum += row[index];
                n++;
            }
        }
        return n > 0 ? sum / n : Double.NaN;
    }

    private double sdAt(double[][] values, int index, double mean) {
        if (!Double.isFinite(mean)) {
            return Double.NaN;
        }
        double sum = 0.0;
        int n = 0;
        for (double[] row : values) {
            if (index < row.length && Double.isFinite(row[index])) {
                double d = row[index] - mean;
                sum += d * d;
                n++;
            }
        }
        return n > 1 ? Math.sqrt(sum / (n - 1)) : 0.0;
    }

    private int countFiniteAt(double[][] values, int index) {
        int n = 0;
        for (double[] row : values) {
            if (index < row.length && Double.isFinite(row[index])) {
                n++;
            }
        }
        return n;
    }

    private double maxFinite(double[] values) {
        double max = Double.NEGATIVE_INFINITY;
        for (double value : values) {
            if (Double.isFinite(value) && value > max) {
                max = value;
            }
        }
        return Double.isFinite(max) ? max : 1.0;
    }

    private double maxFinite(double[][] values, double[] average) {
        double max = Double.NEGATIVE_INFINITY;
        for (double[] row : values) {
            for (double value : row) {
                if (Double.isFinite(value) && value > max) {
                    max = value;
                }
            }
        }
        for (double value : average) {
            if (Double.isFinite(value) && value > max) {
                max = value;
            }
        }
        return Double.isFinite(max) ? max : 1.0;
    }

    private Color sliceColor(int index, int total) {
        if (total <= 1) {
            return Color.BLUE;
        }
        float hue = (float) index / (float) Math.max(1, total);
        return Color.getHSBColor(hue, 0.78f, 0.82f);
    }

    private static class RoiEntry {
        final Roi roi;
        final String name;
        final String entryName;
        String displayName;
        final double centerX;
        final double centerY;
        final double diameterPx;

        RoiEntry(Roi roi, String name, String entryName) {
            this.roi = roi;
            this.name = name;
            this.entryName = entryName;
            double[] center = centerOf(roi);
            this.centerX = center[0];
            this.centerY = center[1];
            Rectangle2D.Double bounds = roi.getFloatBounds();
            this.diameterPx = (bounds.width + bounds.height) / 2.0;
        }

        private static double[] centerOf(Roi roi) {
            Rectangle2D.Double bounds = roi.getFloatBounds();
            double fallbackX = bounds.x + bounds.width / 2.0;
            double fallbackY = bounds.y + bounds.height / 2.0;
            double[] centroid = roi.getContourCentroid();
            if (centroid != null && centroid.length >= 2
                && Double.isFinite(centroid[0]) && Double.isFinite(centroid[1])) {
                return new double[]{centroid[0], centroid[1]};
            }
            return new double[]{fallbackX, fallbackY};
        }
    }

    private static class BatchParameters {
        final RoiEntry roi;
        final String roiZipPath;
        final double objectDiameterMm;
        final double toleranceMm;
        final int maxRadius;
        final int angleSamples;
        final int fftSamples;
        final int startSlice;
        final int endSlice;
        final int sliceStep;
        final boolean confirmPlateaus;
        final boolean showPolarImages;
        final boolean pauseAfterSlice;
        final boolean exportWorkbook;
        PlotWindow[] lastSliceWindows;
        boolean autoRemaining;
        double fixedHighIntersection = Double.NaN;
        double fixedLowIntersection = Double.NaN;

        BatchParameters(RoiEntry roi, String roiZipPath, double objectDiameterMm, double toleranceMm,
                        int maxRadius, int angleSamples, int fftSamples, int startSlice, int endSlice,
                        int sliceStep, boolean confirmPlateaus, boolean showPolarImages,
                        boolean pauseAfterSlice, boolean exportWorkbook) {
            this.roi = roi;
            this.roiZipPath = roiZipPath;
            this.objectDiameterMm = objectDiameterMm;
            this.toleranceMm = toleranceMm;
            this.maxRadius = maxRadius;
            this.angleSamples = angleSamples;
            this.fftSamples = fftSamples;
            this.startSlice = startSlice;
            this.endSlice = endSlice;
            this.sliceStep = sliceStep;
            this.confirmPlateaus = confirmPlateaus;
            this.showPolarImages = showPolarImages;
            this.pauseAfterSlice = pauseAfterSlice;
            this.exportWorkbook = exportWorkbook;
        }
    }

    private static class SliceMeasurement {
        final int slice;
        final SliceResult result;
        final TtfResult ttf;
        final double[] ctf;
        final double contrastHu;

        SliceMeasurement(int slice, SliceResult result, TtfResult ttf, double[] ctf, double contrastHu) {
            this.slice = slice;
            this.result = result;
            this.ttf = ttf;
            this.ctf = ctf;
            this.contrastHu = contrastHu;
        }
    }

    private static class BatchCurves {
        final double[] frequency;
        final double[][] ttfPerSlice;
        final double[][] ctfPerSlice;
        final double[] averageTtf;
        final double[] averageCtf;
        final double[] sdTtf;
        final double[] sdCtf;
        final int[] nTtf;
        final int[] nCtf;
        final int[] sliceNumbers;
        final double nyquist;

        BatchCurves(double[] frequency, double[][] ttfPerSlice, double[][] ctfPerSlice,
                    double[] averageTtf, double[] averageCtf, double[] sdTtf, double[] sdCtf,
                    int[] nTtf, int[] nCtf, int[] sliceNumbers, double nyquist) {
            this.frequency = frequency;
            this.ttfPerSlice = ttfPerSlice;
            this.ctfPerSlice = ctfPerSlice;
            this.averageTtf = averageTtf;
            this.averageCtf = averageCtf;
            this.sdTtf = sdTtf;
            this.sdCtf = sdCtf;
            this.nTtf = nTtf;
            this.nCtf = nCtf;
            this.sliceNumbers = sliceNumbers;
            this.nyquist = nyquist;
        }
    }
}
