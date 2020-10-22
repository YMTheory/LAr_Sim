/*************************************************************************
  @Author: MiaoYu ---> miaoyu@ihep.ac.cn
  @Created Time : Fri Sep 25 15:53:38 2020
  @File Name: LArCollectorHit.hh
 ************************************************************************/

#ifndef LArCollectorHit_h
#define LArCollectorHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
#include "tls.hh"

class LArCollectorHit : public G4VHit
{
    public:
        LArCollectorHit();
        LArCollectorHit(const LArCollectorHit&);
        virtual ~LArCollectorHit();

        // operators
        const LArCollectorHit& operator=(const LArCollectorHit&);
        G4int operator==(const LArCollectorHit&) const;

        inline void* operator new    (size_t);
        inline void operator  delete (void*);

        // methods from base class
        virtual void Draw();
        virtual void Print();

        // Set methods
        void SetTrackID(G4int trackId) { fTrackID = trackId; }

        // Get methods
        G4int  GetTrackID() const;

    private:
        G4double fTrackID;

};


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

typedef G4THitsCollection<LArCollectorHit> LArCollectorHitsCollection;

extern G4ThreadLocal G4Allocator<LArCollectorHit>* LArCollectorHitAllocator;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline void* LArCollectorHit::operator new(size_t)
{
  if(!LArCollectorHitAllocator)
      LArCollectorHitAllocator = new G4Allocator<LArCollectorHit>;
  return (void *) LArCollectorHitAllocator->MallocSingle();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline void LArCollectorHit::operator delete(void *hit)
{
  LArCollectorHitAllocator->FreeSingle((LArCollectorHit*) hit);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

inline G4int LArCollectorHit::GetTrackID() const {
    return fTrackID;
}



#endif
