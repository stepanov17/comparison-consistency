package com.vniim.consistency;

import java.io.File;
import java.io.IOException;
import java.util.concurrent.CountDownLatch;
import javafx.application.Application;
import javafx.application.Platform;
import javafx.embed.swing.SwingFXUtils;
import javafx.scene.Node;
import javafx.scene.Scene;
import javafx.scene.SnapshotParameters;
import javafx.scene.chart.LineChart;
import javafx.scene.chart.NumberAxis;
import javafx.scene.chart.XYChart;
import javafx.scene.control.Tab;
import javafx.scene.control.TabPane;
import javafx.scene.image.WritableImage;
import javafx.scene.layout.BorderPane;
import javafx.stage.Stage;
import javax.imageio.ImageIO;

/**
 * Plot container class.
 */
public class PlotContainer extends Application {

    public static final CountDownLatch latch = new CountDownLatch(2);
    public static PlotContainer plotContainer = null;

    private BorderPane pane;

    private Stage stage;

    /**
     * Get plot container instance.
     * @return the instance of the class.
     */
    public static PlotContainer getInstance() {

        try { latch.await(); }
        catch (InterruptedException e) { e.printStackTrace(); }
        return plotContainer;
    }

    private static void setPlotContainer(PlotContainer c) {

        plotContainer = c;
        latch.countDown();
    }

    /**
     * Constructor.
     */
    public PlotContainer() { setPlotContainer(this); }

    /**
     * Show and save plots to a specified directory.
     * @param plotsData the data ([[x1]; [y1]; [x2]; [y2]; ...)
     * @param outDir the output directory
     */
    public void showAndSavePlots(double plotsData[][], String outDir) {

        Platform.runLater(() -> {

            int nPlots = plotsData.length / 2;
            if (nPlots < 1) { return; }

            TabPane tabPane = new TabPane();

            for (int i = 0; i < nPlots; ++i) {

                Tab tab = new Tab();

                // >>> chart
                double xData[] = plotsData[2 * i];
                double yData[] = plotsData[2 * i + 1];

                String sData;

                double x0, x1, tick, maxX = xData[xData.length - 1];
                if (i < nPlots - 2) {
                    x0 = 1.;
                    x1 = 0.5 * (Math.round(2. * maxX) + 1);
                    tick = 0.5;
                    sData = "k_" + (i + 1);
                } else {
                    x0 = xData[0];
                    x1 = maxX;
                    tick = 0.1 * (x1 - x0);
                    sData = (i == nPlots - 1) ? "uRef" : "xRef";
                }
                tab.setText(sData);

                NumberAxis
                    x = new NumberAxis(x0, x1, tick),
                    y = new NumberAxis(0., 1., 0.1);

                LineChart<Number, Number> plot = new LineChart<>(x, y);
                plot.setAnimated(false);
                plot.setLegendVisible(false);
                XYChart.Series series = new XYChart.Series();

                if (xData.length != yData.length) {
                    throw new IllegalArgumentException("inconsistent data");
                }

                for (int j = 0; j < xData.length; ++j) {
                    series.getData().add(new XYChart.Data(xData[j], yData[j]));
                }
                plot.getData().add(series);
                // <<< chart

                tab.setContent(plot);
                tabPane.getTabs().add(tab);
            }

            pane = new BorderPane();
            pane.setCenter(tabPane);

            Scene scene = new Scene(pane, 700, 500);
            stage.setScene(scene);
            stage.setTitle("k distributions (pdfs)");
            stage.show();


            SnapshotParameters param = new SnapshotParameters();
            for (int i = 0; i < nPlots; ++i) {

                tabPane.getSelectionModel().select(i);

                String sName;
                if (i < nPlots - 2) {
                    sName = "k_" + (i + 1);
                } else  {
                    sName = (i == nPlots - 2) ? "xRef" : "uRef";
                }
                String sFile = outDir + File.separator + sName + ".png";

                Node content = tabPane.getTabs().get(i).getContent();
                //content.requestFocus();
                WritableImage snapShot = content.snapshot(param, null);
                try {
                    File f = new File(sFile);
                    if (f.exists()) {
                        System.out.println("WARNING: overwriting " + sFile);
                    }
                    ImageIO.write(
                        SwingFXUtils.fromFXImage(snapShot, null), "png", f);
                } catch (IOException e) {
                    throw new RuntimeException(e.getMessage());
                }
            }

            tabPane.getSelectionModel().select(0);
        });
    }

    @Override
    public void start(Stage s) {
        stage = s;
        latch.countDown();
    }

    public static void main(String[] args) { Application.launch(args); }
}
