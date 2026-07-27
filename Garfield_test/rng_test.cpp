#include <iomanip>
#include <iostream>

#include "Garfield/Random.hh"
#include "Garfield/RandomEngineRoot.hh"

#include "lgad_tools.hh"

void PrintGlobalSequence(const unsigned long seed) {
  SeedGarfieldRandom(seed);

  for (int i = 0; i < 5; ++i) {
    std::cout << std::setprecision(17)
              << Garfield::RndmUniform() << " ";
  }
  std::cout << "\n";
}

void PrintDirectSequence(const unsigned long seed) {
  Garfield::RandomEngineRoot engine;
  engine.SetSeed(static_cast<unsigned int>(seed));

  for (int i = 0; i < 5; ++i) {
    std::cout << std::setprecision(17)
              << engine.Draw() << " ";
  }
  std::cout << "\n";
}

int main() {
  std::cout << "Direct engine, seed 12345:\n";
  PrintDirectSequence(12345);
  PrintDirectSequence(12345);

  std::cout << "\nInstalled Garfield engine, seed 12345:\n";
  PrintGlobalSequence(12345);
  PrintGlobalSequence(12345);

  std::cout << "\nInstalled Garfield engine, seed 12346:\n";
  PrintGlobalSequence(12346);
  

  return 0;
}