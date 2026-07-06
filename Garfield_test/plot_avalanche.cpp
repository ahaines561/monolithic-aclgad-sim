#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>

#include <TCanvas.h>
#include <TAxis.h>
#include <TH1F.h>
#include <TApplication.h>

void plot_avalanche(int argc, char** argv) {
    const std::string file_path = argc > 1 ? argv[1]
        : "/home/ahaines561/HEP/MAS/Garfield_test/build/eh_sizes.txt";

    std::ifstream file(file_path);
    std::vector<double> sizes;
    double size;

    while (file >> size) {
        sizes.push_back(size);
    }
    file.close();

    double max_size = *std::max_element(sizes.begin(), sizes.end());
    TCanvas* c1 = new TCanvas("canvas", "Avalanche size distribution", 800, 600);
    TH1F* hist = new TH1F("hist", "Avalanche size distribution; avalanche size (electrons per injected electron);Events", 50, 0, max_size * 1.1);
    for (double s : sizes) {
        hist->Fill(s);
    }
    hist->SetLineColor(kBlue);
    hist->SetFillColorAlpha(kBlue, 0.3);
    hist->SetLineWidth(2);
    hist->Draw("HIST");
    c1->Update();
    c1->Modified();
    c1->SaveAs("avalanche_histogram.png");
}

int main(int argc, char** argv) {
    TApplication app("app", &argc, argv);
    
    plot_avalanche(argc, argv);
    
    app.Run();
    return 0;
}