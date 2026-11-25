#pragma once

#include "Clock.h"
#include "Host/WithinHost/WHInterface.h"
#include "PkPd/LSTMModel.h"

using namespace std;

namespace OpenMalaria {
namespace UnitTest {

/**
 * State behind the WHInterface Kind::Mock backend. The mock class
 * has no actual infections: just set totalDensity as you please; nTreatments
 * counts treatment()/treatSimple()/treatPkPd() calls; lastTimeLiver/Blood are
 * the last time-spans used by the simple treatment model (never() if not
 * used); pkpd lists medications and drugs in the body.
 */
struct WHMockState {
    double totalDensity = numeric_limits<double>::quiet_NaN();
    int nTreatments = 0;
    SimTime lastTimeLiver = never(), lastTimeBlood = never();
    LSTMModel pkpd;
};

/**
 * Create a Kind::Mock within-host model for a human (see
 * UnittestUtil::setHumanWH). Each mock owns its independent state.
 */
unique_ptr<WHInterface> makeWHMock(RunContext& ctx);
WHMockState& whMockState(WHInterface& mock);

}
}
