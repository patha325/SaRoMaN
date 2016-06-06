//
//  MindDigiHit.hpp
//  
//
//  Created by Ryan Bayes on 2016-05-31.
//
//

#ifndef MindDigiHit_hpp
#define MindDigiHit_hpp

#include <stdio.h>

#endif /* MindDigiHit_hpp */

class MindDigiHit: public G4VDigi
{
public:
    /// Constructor
    MindDigiHit() {}
    /// Destuctor
    ~MindDigiHit() {}
    /// Copy-constructor
    MindDigiHit(const MindDigiHit&);
    /// Assignment operator
    const MindDigiHit& operator=(const MindDigiHit&);
    /// Equality operator
    G4int operator==(const MindDigiHit&) const;
    
    /// Memory allocation
    inline void* operator new(size_t);
    inline void  operator delete(void*);
    
public:
    /// Setters and Getters
    ///
    G4int GetTrackID() { return _track_id; }
    void  SetTrackID(G4int id) { _track_id = id; }
    ///
    G4ThreeVector GetPosition() { return _position; }
    void SetPosition(G4ThreeVector p) { _position = p; }
    ///
    G4double GetEnergyDeposit() { return _energy_dep; }
    void SetEnergyDeposit(G4double e) { _energy_dep = e; }
    //
    G4double GetHitTime() { return _time; }
    void SetHitTime(G4double t) { _time = t; }
    
    G4ThreeVector GetBarTranslation() { return _barnumber; }
    void SetBarTranslation(G4ThreeVector b) { _barnumber = b; }
    
    G4String GetModule() { return _module; }
    void SetModule(G4String m) { _module = m; }
    
    G4int GetBarOrientation() { return _isYbar; }
    void SetBarOrientation(G4int i) { _isYbar = i; }
    
    G4int GetBarNumber() { return _barNumber; }
    void SetBarNumber(G4int i) { _barNumber = i; }
    
    G4int GetChannelID() { return _channelid; }
    void SetChannelID(G4int i) { _channelid = i; }
    
    G4int GetBoardNumber() { return _boardnumber; }
    void SetBoardNumber(G4int i) { _boardnumber; }
    
    G4double GetTimeOverThreshold() { return _timeoverth; }
    void SetTimeOverThreshold(G4double a) { _timeoverth = a; }
    
private:
    G4int _track_id;
    G4double _energy_dep;
    G4ThreeVector _position;
    G4ThreeVector _barnumber;
    G4double _time;
    G4String _module;
    G4int _channelid;
    G4int _boardnumber;
    G4double _timeoverth;
};