#ifndef Hmod_WithinHost_Mock
#define Hmod_WithinHost_Mock

#include "Global.h"
#include "Host/WithinHost/WHInterface.h"
#include "PkPd/LSTMModel.h"

using namespace std;

namespace OM {
namespace UnitTest {

/**
 * State behind the WithinHost::WHInterface Kind::Mock backend. The mock class
 * has no actual infections: just set totalDensity as you please; nTreatments
 * counts treatment()/treatSimple()/treatPkPd() calls; lastTimeLiver/Blood are
 * the last time-spans used by the simple treatment model (sim::never() if not
 * used); pkpd lists medications and drugs in the body.
 */
struct WHMockState {
    double totalDensity = numeric_limits<double>::quiet_NaN();
    int nTreatments = 0;
    SimTime lastTimeLiver = sim::never(), lastTimeBlood = sim::never();
    PkPd::LSTMModel pkpd;
};

/**
 * Create a Kind::Mock within-host model for a human (see
 * UnittestUtil::setHumanWH), installing fresh state and the hook table on
 * ctx.ext. Returns the model; the state is readable via
 * ctx.ext->whMock->state.
 */
unique_ptr<WithinHost::WHInterface> makeWHMock(RunContext& ctx);

}
}
#endif
