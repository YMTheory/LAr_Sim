#ifndef LArActionInitialization_h
#define LArActionInitialization_h 1

#include "G4VUserActionInitialization.hh"

/// Action initialization class.

class LArActionInitialization : public G4VUserActionInitialization
{
  public:
    LArActionInitialization();
    virtual ~LArActionInitialization();

    virtual void Build() const;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

    
