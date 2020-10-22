/*************************************************************************
 @Author: MiaoYu ---> miaoyu@ihep.ac.cn
 @Created Time : Fri Sep 25 15:58:10 2020
 @File Name: LArCollectorHit.cc
 ************************************************************************/

#include "LArCollectorHit.hh"

G4ThreadLocal G4Allocator<LArCollectorHit>* LArCollectorHitAllocator = 0;

LArCollectorHit::LArCollectorHit()
    : G4VHit(),
      fTrackID(0.)
{}

LArCollectorHit::~LArCollectorHit()
{;}


LArCollectorHit::LArCollectorHit( const LArCollectorHit& right)
    : G4VHit()
{
    fTrackID = right.fTrackID;
}


const LArCollectorHit& LArCollectorHit::operator=(const LArCollectorHit& right)
{
    fTrackID = right.fTrackID;

    return *this;
}

G4int LArCollectorHit::operator==(const LArCollectorHit& right) const
{
    return (  this == &right ) ? 1 : 0;
}
void LArCollectorHit::Draw()
{

}

void LArCollectorHit::Print()
{
    G4cout << " trackID: " << fTrackID << G4endl;
}





