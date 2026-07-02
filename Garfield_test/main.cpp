#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <sstream>
#include <fstream>
#include <cmath>
#include <limits>
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

// read mesh bounds
bool GetMeshBounds(const std::string& path, double& xmin, double& ymin,
                   double& xmax, double& ymax) {
    std::ifstream fh(path);
    if (!fh) return false;
    xmin = ymin = std::numeric_limits<double>::max();
    xmax = ymax = std::numeric_limits<double>::lowest();
    std::string line;
    bool found = false;
    while (std::getline(fh, line)) {
        if (line.empty() || line[0] != 'c') continue;
        std::istringstream ss(line);
        std::string tag; double id, x, y;
        if (!(ss >> tag >> id >> x >> y) || tag != "c") continue;
        xmin = std::min(xmin, x); xmax = std::max(xmax, x);
        ymin = std::min(ymin, y); ymax = std::max(ymax, y);
        found = true;
    }
    if (!found) return false;
    xmin *= 1.e-4; xmax *= 1.e-4; ymin *= 1.e-4; ymax *= 1.e-4;
    return true;
}

// read region->material map from the 'r id material' records
std::map<int, std::string> GetRegionMaterials(const std::string& path) {
    std::map<int, std::string> mats;
    std::ifstream fh(path);
    if (!fh) return mats;
    std::string line;
    while (std::getline(fh, line)) {
        if (line.empty() || line[0] != 'r') continue;
        std::istringstream ss(line);
        std::string tag; int id; std::string material;
        if (!(ss >> tag >> id >> material) || tag != "r") continue;
        mats[id] = material;
    }
    return mats;
}

int main(int argc, char* argv[]) {
    TApplication app("app", &argc, argv);

    const std::string file =
        argc > 1 ? argv[1]
                 : "/home/ahaines561/HEP/MAS/Silvaco_dat/lgad.sta";

    // signal timing granularity
    const int nBins = 1000;
    const double tmax = 5.0;
    const double dt = tmax / nBins;
    const int nEvents = 10000;
    int totElectronsEvent0 = 0;
    int totHolesEvent0 = 0;
    int nPrimaryEvent0 = 0;
    int nPrimaryE = 0;
    int eDone = 0;
    bool signalDone = false;
    double avgGain = 0.0;
    double t = 0.0;
    double current = 0.0;

    const bool readoutAtTop = true;
    const double tauAC = 0.5;

    MediumSilicon si;
    si.SetTemperature(293.15);
    si.SetImpactIonisationModelOkutoCrowell();
    // si.SetImpactIonisationModelMassey();
    si.SetSaturationVelocityModelCanali();

    ComponentTcad2d field;
    if (!field.InitialiseSilvaco(file)) {
        std::cerr << "InitialiseSilvaco FAILED\n";
        return 1;
    }

    // read geometry
    double xmin, ymin, xmax, ymax;
    if (!GetMeshBounds(file, xmin, ymin, xmax, ymax)) {
        std::cerr << "GetMeshBounds FAILED\n";
        return 1;
    }
    const auto mats = GetRegionMaterials(file);

    // assign si to mat
    int nSiRegions = 0;
    const std::size_t nReg = field.GetNumberOfRegions();
    for (std::size_t i = 0; i < nReg; ++i) {
        std::string material = mats.count(i) ? mats.at(i) : "";
        bool isSilicon = (material == "1" || material == "3");
        if (isSilicon) { field.SetMedium(i, &si); ++nSiRegions; }
    }

    // sanity check
    std::cout << "\n device summary \n";
    std::cout << "file   : " << file << "\n";
    std::cout << "bounds : x[" << xmin << ", " << xmax << "]  y[" << ymin
              << ", " << ymax << "] cm\n";
    std::cout << "regions: " << nReg << " total, " << nSiRegions << " silicon\n";
    for (std::size_t i = 0; i < nReg; ++i) {
        std::string material = mats.count(i) ? mats.at(i) : "";
        std::cout << "  region " << i << " material='" << material << "'"
                  << ((material == "1" || material == "3") ? " -> silicon" : "")
                  << "\n";
    }
    if (nSiRegions == 0) {
        std::cerr << "WARNING: no silicon regions found - check 'r' records!\n";
    }

    // -geometry from bounds
    const double ySi_Top = ymin;
    const double ySi_Bot = ymax;
    const double xCenter = 0.5 * (xmin + xmax);
    const double thickness = ySi_Bot - ySi_Top;

    ComponentUser wpad;
    wpad.SetArea(xmin, ymin, -1., xmax, ymax, 1.);
    wpad.SetMedium(&si);
    // linear weighting potential; readout at ySi_Top (flip if readoutAtTop=false)
    wpad.SetWeightingPotential(
        [ySi_Top, ySi_Bot, readoutAtTop](const double, const double y, const double) {
            const double yc = std::min(std::max(y, ySi_Top), ySi_Bot);
            const double w = (ySi_Bot - yc) / (ySi_Bot - ySi_Top);
            return readoutAtTop ? w : (1.0 - w);
        }, "pad");
    wpad.SetWeightingField(
        [ySi_Top, ySi_Bot, readoutAtTop](const double, const double y, const double,
                     double& wx, double& wy, double& wz) {
            wx = 0.;
            const double mag = 1.0 / (ySi_Bot - ySi_Top);
            wy = (y >= ySi_Top && y <= ySi_Bot) ? (readoutAtTop ? mag : -mag) : 0.;
            wz = 0.;
        }, "pad");

    // Sensor
    Sensor sensor;
    sensor.AddComponent(&field);
    sensor.AddElectrode(&wpad, "pad");
    sensor.SetTimeWindow(0., dt, nBins);
    sensor.SetArea(xmin, ymin, -1., xmax, ymax, 1.);

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

    // start the track near the back face, drifting down through the depth
    const double yStart = ySi_Bot - 0.1 * thickness;

    std::cout << "Shooting " << nEvents << " pions...\n";

    for (int i = 0; i < nEvents; ++i) {
        track.NewTrack(xCenter, yStart, 0., 0., 0., -1., 0.);

        const auto& clusters = track.GetClusters();
        nPrimaryE = 0;
        for (const auto& cluster : clusters) {
            nPrimaryE += cluster.electrons.size();
        }
        hLandau->Fill(nPrimaryE);
        // generate the waveform on the first event with enough primaries
        if (!signalDone && nPrimaryE > 1000){
            signalDone = true;
            std::cout << "\nRunning full Avalanche on event " << i << " (primaries: " << nPrimaryE << ")...\n";
            nPrimaryEvent0 = nPrimaryE;
            eDone = 0;
            for (const auto& cluster : clusters){
                for (const auto& e : cluster.electrons) {
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
            std::cout << "\nEvent avalanche complete. Landau generation continues ...\n";
        }
        if (i > 0 && i % 500 == 0) {
            std::cout << "  Processed " << i << " total events...\r" << std::flush;
        }
    }
    std::cout << "\nFinished " << nEvents << " events.\n";
    if (nPrimaryEvent0 > 0) {
        avgGain = static_cast<double>(totElectronsEvent0) / nPrimaryEvent0;
    }
    std::cout << "Gain: " << avgGain << "\n";

    // save the raw (DC) signal before AC coupling
    std::ofstream rawfile("dc_signal.txt");
    if (rawfile.is_open()){
        rawfile << "Time_ns\tCurrent\n";
        for (int i = 0; i < nBins; ++i) {
            t = (i + 0.5) * dt;
            rawfile << t << "\t" << sensor.GetSignal("pad", i) << "\n";
        }
        rawfile.close();
    }

    // plot the raw DC pulse
    TCanvas* c0 = new TCanvas("c0", "DC Pulse (raw)", 800, 600);
    ViewSignal svDC(&sensor);
    svDC.SetCanvas(c0);
    svDC.PlotSignal("pad", "t");
    c0->Update();
    c0->SaveAs("dc_pulse.png");

    // AC coupling: bipolar CR high-pass kernel -> bipolar AC-LGAD shape
    sensor.SetTransferFunction(
        [tauAC](double tt) {
            return (tt < 0.) ? 0.
                             : (1.0 - tt / tauAC) * std::exp(-tt / tauAC) / tauAC;
        });
    sensor.ConvoluteSignals();

    // stats + AC signal
    std::ofstream outfile("avalanche_results.txt");
    if (outfile.is_open()){
        outfile << "--- Event 0 AC-LGAD Summary ---\n";
        outfile << "Primary pairs: " << nPrimaryEvent0 << "\n";
        outfile << "Total electrons: " << totElectronsEvent0 << "\n";
        outfile << "Total holes: " << totHolesEvent0 << "\n";
        outfile << "Gain: " << avgGain << "\n";
        outfile << "\n--- AC Signal Data ---\nTime_ns\tCurrent\n";
        for (int i = 0; i < nBins; ++i) {
            t = (i + 0.5) * dt;
            current = sensor.GetSignal("pad", i);
            outfile << t << "\t" << current << "\n";
        }
        outfile.close();
    }

    // plot the AC-coupled (bipolar) pulse
    TCanvas* c1 = new TCanvas("c1", "AC-LGAD Pulse", 800, 600);
    ViewSignal sv(&sensor);
    sv.SetCanvas(c1);
    sv.PlotSignal("pad", "t");
    c1->Update();
    c1->SaveAs("ac_pulse.png");

    TCanvas* c2 = new TCanvas("c2", "Landau Distribution", 800, 600);
    hLandau->Draw();
    c2->Update();
    c2->SaveAs("landau_histogram.png");

    app.Run();
    return 0;
}