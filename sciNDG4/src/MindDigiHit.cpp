//
//  MindDigiHit.cpp
//  
//
//  Created by Ryan Bayes on 2016-05-31.
//
//

#include "MindDigiHit.h"
#include <G4Allocator.hh>

G4Allocator<MindDigiHit> MindDigiHitAllocator;

MindDigiHit::MindDigiHit(const MindDigiHit& right): G4VDigi()
{
    _track_id   = right._track_id;
    _energy_dep = right._energy_dep;
    _position   = right._position;
    _time       = right._time;
    _barnumber  = right._barnumber;
    _module     = right._module;
    _isYbar     = right._isYbar;
    // Raw information expected to be returned by the digitization
    _boardnumber  = right._boardnumber;
    _channelid    = right._channelid;
    _timeoverth   = right._timeoverth;
}

const MindDigiHit& MindDigiHit::operator=(const MindDigiHit& right)
{
    _track_id   = right._track_id;
    _energy_dep = right._energy_dep;
    _position   = right._position;
    _time       = right._time;
    _barnumber  = right._barnumber;
    _module     = right._module;
    _isYbar     = right._isYbar;
    // Raw information expected to be returned by the digitization
    _boardnumber  = right._boardnumber;
    _channelid    = right._channelid;
    _timeoverth   = right._timeoverth;
}

G4int MindDigiHit::operator==(const MindDigiHit& right) const
{
    return (this==&right) ? 1 : 0;
}

std::bitset<32> MindDigiHit::