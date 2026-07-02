#include <iostream>
#include <string>
#include <map>
#include <sstream>
#include <fstream>
#include <cmath>
#include <limits>

#include "Garfield/MediumSilicon.hh"
#include "Garfield/ComponentTcad2d.hh"

using namespace Garfield;

bool GetMeshBounds(const std::string& path, double& xmin, double& ymin,
                   double& xmax, double& ymax) {
    std::ifstream fh(path);
    if (!fh) return false;
    xmin = ymin = std::numeric_limits<double>::max();
    xmax = ymax = std::numeric_limits<double>::lowest();
    std::string line; bool found = false;
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
std::map<int, std::string> GetRegionMaterials(const std::string& path) {
    std::map<int, std::string> mats;
    std::ifstream fh(path); if (!fh) return mats;
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
    const std::string file =
        argc > 1 ? argv[1] : "/home/ahaines561/HEP/MAS/Silvaco_dat/lgad.sta";
    MediumSilicon si; si.SetTemperature(293.15);
    ComponentTcad2d field;
    if (!field.InitialiseSilvaco(file)) { std::cerr << "load failed\n"; return 1; }
    double xmin, ymin, xmax, ymax;
    GetMeshBounds(file, xmin, ymin, xmax, ymax);
    const auto mats = GetRegionMaterials(file);
    const std::size_t nReg = field.GetNumberOfRegions();
    for (std::size_t i = 0; i < nReg; ++i) {
        std::string m = mats.count(i) ? mats.at(i) : "";
        if (m == "1" || m == "3") field.SetMedium(i, &si);
    }

    // fine grid over the whole device; track max |E| in silicon (medium != null)
    double bestE = 0, bestX = 0, bestY = 0;
    const int nx = 200, ny = 200;
    for (int ix = 0; ix <= nx; ++ix) {
        for (int iy = 0; iy <= ny; ++iy) {
            double x = xmin + (xmax - xmin) * ix / nx;
            double y = ymin + (ymax - ymin) * iy / ny;
            double ex, ey, ez, v; Medium* m = nullptr; int st = 0;
            field.ElectricField(x, y, 0., ex, ey, ez, v, m, st);
            if (!m) continue;  // only silicon (has a medium)
            double emag = std::sqrt(ex*ex + ey*ey + ez*ez);
            if (emag > bestE) { bestE = emag; bestX = x; bestY = y; }
        }
    }
    std::cout << "Max |E| in silicon: " << bestE << " V/cm at x = "
              << bestX*1e4 << " um, y = " << bestY*1e4 << " um\n";

    // also show the vertical field profile down that x
    std::cout << "\nField down x = " << bestX*1e4 << " um:\n y[um]\t|E|[V/cm]\n";
    for (int iy = 0; iy <= 60; ++iy) {
        double y = ymin + (ymax - ymin) * iy / 60.0;
        double ex, ey, ez, v; Medium* m = nullptr; int st = 0;
        field.ElectricField(bestX, y, 0., ex, ey, ez, v, m, st);
        double emag = std::sqrt(ex*ex + ey*ey + ez*ez);
        std::cout << y*1e4 << "\t" << emag << (m ? "" : "  (no medium)") << "\n";
    }
    return 0;
}