#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <fstream>
#include <TApplication.h>
#include <TCanvas.h>
#include <TH1F.h>

#include "Garfield/MediumSilicon.hh"
#include "Garfield/ComponentTcad2d.hh"
#include "Garfield/ComponentUser.hh"
#include "Garfield/Sensor.hh"
#include "Garfield/TrackHeed.hh"
#include "Garfield/AvalancheMC.hh"
#include "Garfield/ViewSignal.hh"

using namespace Garfield;

int main(int argc, char* argv[]) {
    TApplication app("app", &argc, argv);

    const std::string file =
        argc > 1 ? argv[1]
                 : "/home/ahaines561/HEP/MAS/Silvaco_dat/"
                   "30kev5umsp_5um_X26.sta";

    constexpr double yBot = -3.60e-4;
    constexpr double ySi_Top = 0.0e-4;
    constexpr double ySi_Bot = 50.50e-4;
    const int nBins = 1000;
    const double tmax = 5.0;
    const double dt = tmax / nBins;
    const int nEvents = 10000;
    int totElectronsEvent0 = 0;
    int totHolesEvent0 = 0;
    int nPrimaryEvent0 = 0;
    int nPrimaryE = 0;
    int eDone = 0;
    double avgGain = 0.0;
    double t = 0.0;
    double current = 0.0;

    MediumSilicon si;
    si.SetTemperature(293.15);
    // si.SetImpactIonisationModelOkutoCrowell();
    si.SetImpactIonisationModelMassey();
    si.SetSaturationVelocityModelCanali();

    ComponentTcad2d field;
    if (!field.InitialiseSilvaco(file)) {
        std::cerr << "InitialiseSilvaco FAILED\n";
        return 1;
    }
    
    for (int i = 1; i <= 13; ++i) {
        field.SetMedium(i, &si);
    }
    ComponentUser wpad;
    // SetArea: (xMin, yMin, zMin, xMax, yMax, zMax)
    wpad.SetArea(0., yBot, -1., 250.e-4, ySi_Bot, 1.);
    wpad.SetMedium(&si);
    //weighting potential
    wpad.SetWeightingPotential(
        [ySi_Top, ySi_Bot](const double, const double y, const double) {
            const double yc = std::min(std::max(y, ySi_Top), ySi_Bot);
            return (ySi_Bot - yc) / (ySi_Bot - ySi_Top);
        }, "pad");
    wpad.SetWeightingField(
        [ySi_Top, ySi_Bot](const double, const double y, const double,
                     double& wx, double& wy, double& wz) {
            wx = 0.;
            wy = (y >= ySi_Top && y <= ySi_Bot) ? 1.0 / (ySi_Bot - ySi_Top) : 0.;
            wz = 0.;
        }, "pad");
    // Sensor
    Sensor sensor;
    sensor.AddComponent(&field);
    sensor.AddElectrode(&wpad, "pad");
    sensor.SetTimeWindow(0., dt, nBins);
    sensor.SetArea(0., ySi_Top, -1., 250.e-4, ySi_Bot, 1.);

    // MIP
    TrackHeed track(&sensor);
    track.SetParticle("pion");
    track.SetMomentum(180.e9);
    track.EnableDeltaElectronTransport();

    // Avalanche
    AvalancheMC aval(&sensor);
    aval.SetTimeSteps(0.001); 
    aval.EnableSignalCalculation();
    aval.EnableAvalancheSizeLimit(100); 
    aval.SetTimeWindow(0., 10.0); 

    TH1F* hLandau = new TH1F("hLandau", "Primary Charge Distribution;Total Primary Electrons;Frequency", 100, -100, 5000);

    // run the Simulation
    std::cout << "Shooting " << nEvents << " pions...\n";

    for (int i = 0; i < nEvents; ++i) {
        // Shoot the track straight up the middle of the 250um wide sensor
        track.NewTrack(125.e-4, 45.e-4, 0., 1.0, 0., -1., 0.);

        const auto& clusters = track.GetClusters();
        nPrimaryE = 0;
        for (const auto& cluster : clusters) {
            nPrimaryE += cluster.electrons.size();
        }
        hLandau->Fill(nPrimaryE);
        if (i == 0){
            std::cout << "\nRunning full Avalanche on Event 0 to generate waveform...\n";
            nPrimaryEvent0 = nPrimaryE;
            eDone = 0;
            // transport primaries
            for (const auto& cluster : clusters){
                for (const auto& e : cluster.electrons) {
                    // avalanche the electron + transport the primary and avalanche holes
                    aval.AvalancheElectronHole(e.x, e.y, e.z, e.t);
                    totElectronsEvent0 += aval.GetElectrons().size();
                    totHolesEvent0 += aval.GetHoles().size();
                    eDone++;
                    if (eDone % 50 == 0 || eDone == nPrimaryE) {
                        std::cout << "  Avalanched " << eDone << "/" << nPrimaryE 
                                  << " primary electrons (" 
                                  << static_cast<int>(100.0 * eDone / nPrimaryE) << "%)\r" 
                                  << std::flush;
                    }
                }
            }
            std::cout << "\nEvent 0 Avalanche complete. Landau generation next ...\n";
        }
        if (i > 0 && i % 500 == 0) {
            std::cout << "  Processed " << i << " total events...\r" << std::flush;
        }
    }
    std::cout << "\nFinished 10,000 events.\n";
    if (nPrimaryEvent0 > 0) {
        avgGain = static_cast<double>(totElectronsEvent0) / nPrimaryEvent0;
    }
    //avalanche + landau stats
    std::ofstream outfile("avalanche_results.txt");
    if (outfile.is_open()){
        outfile << "--- Event 0 AC-LGAD Summary ---\n";
        outfile << "Primary pairs: " << nPrimaryEvent0 << "\n";
        outfile << "Total electrons: " << totElectronsEvent0 << "\n";
        outfile << "Total holes: " << totHolesEvent0 << "\n";
        outfile << "Gain: " << avgGain << "\n";
        outfile << "\n--- Signal Data ---\nTime_ns\tCurrent\n";
        
        for (int i = 0; i < nBins; ++i) {
            t = (i + 0.5) * dt; 
            current = sensor.GetSignal("pad", i);
            outfile << t << "\t" << current << "\n";
        }
        outfile.close();
    }

    // Plotting
    TCanvas* c1 = new TCanvas("c1", "ACLGAD Pulse Signal", 800, 600);
    ViewSignal sv(&sensor);
    sv.SetCanvas(c1);
    sv.PlotSignal("pad", "t");
    c1->Update();
    c1->SaveAs("lgad_raw_signal.png");

    TCanvas* c2 = new TCanvas("c2", "Landau Distribution", 800, 600);
    hLandau->Draw();
    c2->Update();
    c2->SaveAs("landau_histogram.png");

    app.Run();
    return 0;
}