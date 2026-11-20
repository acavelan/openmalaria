#include "Host/WithinHost/WHInterface.h"
#include "util/checkpoint.h"
#include "core/RunExtensions.h"
#include "Host/WithinHost/Pathogenesis/PathogenesisModel.h"
#include "Host/WithinHost/Treatments.h"
#include "util/ModelOptions.h"
#include "util/errors.h"

#include <cmath>

namespace OM {
namespace WithinHost {

using namespace OM::util;

bool reportInfectionsByGenotype = false;

bool opt_vivax_simple = false,
        opt_dummy_whm = false, opt_empirical_whm = false,
        opt_molineaux_whm = false, opt_penny_whm = false,
        opt_common_whm = false;

// -----  static functions  -----


TreatmentId WHInterface::addTreatment(RunContext& ctx, const scnXml::TreatmentOption& desc){
    return Treatments::addTreatment( ctx, desc );
}

unique_ptr<WHInterface> WHInterface::createWithinHostModel(RunContext& ctx, LocalRng& rng, double comorbidityFactor) {
    if( opt_vivax_simple ) {
        return unique_ptr<WHInterface>(new WHInterface( ctx, rng, comorbidityFactor, Kind::Vivax ));
    } else if( opt_common_whm ) {
        return unique_ptr<WHInterface>(new WHInterface( ctx, rng, comorbidityFactor, Kind::Common ));
    } else {
        return unique_ptr<WHInterface>(new WHInterface( ctx, rng, comorbidityFactor, Kind::Descriptive ));
    }
}


// -----  Constructors  -----

// parity: order — each ctor* fragment draws exactly where the matching old
// class constructor did (falciparum immunity block before the per-kind body).
WHInterface::WHInterface(RunContext& ctx, LocalRng& rng, double comorbidityFactor, Kind kind) :
    ctx(ctx), m_kind(kind)
{
    switch( m_kind ){
    case Kind::Vivax: ctorVivax( rng, comorbidityFactor ); break;
    case Kind::Common: ctorCommon( rng, comorbidityFactor ); break;
    case Kind::Descriptive: ctorDescriptive( rng, comorbidityFactor ); break;
    case Kind::Mock: break;    // unittest backend: no state, no draws
    }
}

WHInterface::~WHInterface(){
    switch( m_kind ){
    case Kind::Vivax: dtorVivax(); break;
    case Kind::Common: dtorCommon(); break;
    case Kind::Descriptive: dtorDescriptive(); break;
    case Kind::Mock: break;
    }
}

WHInterface::WHInterface(WHInterface&&) = default;


// -----  Kind dispatch  -----
//
// The Kind::Mock branch calls through the hook table on ctx.ext (installed by
// unittest/WHMock.cpp); other Mock branches reproduce the old WHMock's
// unimplemented_exception() behaviour directly. Production never constructs
// the Mock kind.


double WHInterface::probTransmissionToMosquito(vector<double> &probTransGenotype_i, vector<double> &probTransGenotype_l) const{
    switch( m_kind ){
    case Kind::Vivax: return probTransmissionToMosquito_vivax( probTransGenotype_i, probTransGenotype_l );
    case Kind::Mock: throw util::unimplemented_exception( "not needed in unit test" );
    default: return probTransmissionToMosquito_falciparum( probTransGenotype_i, probTransGenotype_l );
    }
}

bool WHInterface::summarize(Host::Human& human) const{
    switch( m_kind ){
    case Kind::Vivax: return summarize_vivax( human );
    case Kind::Common: return summarize_common( human );
    case Kind::Descriptive: return summarize_descriptive( human );
    case Kind::Mock: throw util::unimplemented_exception( "not needed in unit test" );
    }
    return false;   // unreachable
}

void WHInterface::importInfection(LocalRng& rng){
    switch( m_kind ){
    case Kind::Vivax: importInfection_vivax( rng ); break;
    case Kind::Common: importInfection_common( rng ); break;
    case Kind::Descriptive: importInfection_descriptive( rng ); break;
    case Kind::Mock: throw util::unimplemented_exception( "not needed in unit test" );
    }
}

void WHInterface::treatment( Host::Human& human, TreatmentId treatId ){
    if( m_kind == Kind::Mock ){
        const RunExtensions::WHMockHooks* hooks = ctx.ext->whMock;
        hooks->treatment( hooks->state );
        return;
    }
    treatment_falciparum( human, treatId );    // Vivax kind uses the same code
}

void WHInterface::optionalPqTreatment( Host::Human& human ){
    if( m_kind == Kind::Vivax ) return optionalPqTreatment_vivax( human );
    if( m_kind == Kind::Mock ) throw util::unimplemented_exception( "not needed in unit test" );
    // No PQ treatment for falciparum in current models.
}

bool WHInterface::treatSimple( Host::Human& human, SimTime timeLiver, SimTime timeBlood ){
    switch( m_kind ){
    case Kind::Vivax: return treatSimple_vivax( human, timeLiver, timeBlood );
    case Kind::Mock: {
        const RunExtensions::WHMockHooks* hooks = ctx.ext->whMock;
        return hooks->treatSimple( hooks->state, timeLiver, timeBlood );
    }
    default: return treatSimple_falciparum( human, timeLiver, timeBlood );
    }
}

void WHInterface::treatPkPd(size_t schedule, size_t dosages, double age, double delay_d){
    switch( m_kind ){
    case Kind::Vivax: treatPkPd_vivax( schedule, dosages, age, delay_d ); break;
    case Kind::Common: treatPkPd_common( schedule, dosages, age, delay_d ); break;
    case Kind::Descriptive: treatPkPd_descriptive( schedule, dosages, age, delay_d ); break;
    case Kind::Mock: {
        const RunExtensions::WHMockHooks* hooks = ctx.ext->whMock;
        hooks->treatPkPd( hooks->state, ctx, schedule, dosages, age, delay_d );
    }
    }
}

void WHInterface::update(Host::Human &human, LocalRng& rng, int &nNewInfs_i, int &nNewInfs_l, 
        vector<double>& genotype_weights_i, vector<double>& genotype_weights_l, double ageInYears)
{
    switch( m_kind ){
    case Kind::Vivax: update_vivax( human, rng, nNewInfs_i, nNewInfs_l, genotype_weights_i, genotype_weights_l, ageInYears ); break;
    case Kind::Common: update_common( human, rng, nNewInfs_i, nNewInfs_l, genotype_weights_i, genotype_weights_l, ageInYears ); break;
    case Kind::Descriptive: update_descriptive( human, rng, nNewInfs_i, nNewInfs_l, genotype_weights_i, genotype_weights_l, ageInYears ); break;
    case Kind::Mock: throw util::unimplemented_exception( "not needed in unit test" );
    }
}

double WHInterface::getTotalDensity() const{
    if( m_kind == Kind::Vivax ) return getTotalDensity_vivax();
    if( m_kind == Kind::Mock ){
        const RunExtensions::WHMockHooks* hooks = ctx.ext->whMock;
        return hooks->getTotalDensity( hooks->state );
    }
    return totalDensity;
}

bool WHInterface::diagnosticResult( LocalRng& rng, const Diagnostic& diagnostic ) const{
    if( m_kind == Kind::Vivax ) return diagnosticResult_vivax( rng, diagnostic );
    if( m_kind == Kind::Mock ){
        const RunExtensions::WHMockHooks* hooks = ctx.ext->whMock;
        return hooks->diagnosticResult( hooks->state, rng, diagnostic );
    }
    return diagnosticResult_falciparum( rng, diagnostic );
}

Pathogenesis::StatePair WHInterface::determineMorbidity( Host::Human& human, double ageYears, bool isDoomed ){
    if( m_kind == Kind::Vivax ) return determineMorbidity_vivax( human, ageYears, isDoomed );
    if( m_kind == Kind::Mock ) throw util::unimplemented_exception( "not needed in unit test" );
    return determineMorbidity_falciparum( human, ageYears, isDoomed );
}

void WHInterface::clearImmunity(){
    if( m_kind == Kind::Vivax ) return clearImmunity_vivax();
    if( m_kind == Kind::Mock ) throw util::unimplemented_exception( "not needed in unit test" );
    clearImmunity_falciparum();
}

double WHInterface::getCumulative_h() const{
    if( m_kind == Kind::Vivax ) return getCumulative_h_vivax();
    if( m_kind == Kind::Mock ) throw util::unimplemented_exception( "not needed in unit test" );
    return imm_h_ref();
}

double WHInterface::getCumulative_Y() const{
    if( m_kind == Kind::Vivax ) return getCumulative_Y_vivax();
    if( m_kind == Kind::Mock ) throw util::unimplemented_exception( "not needed in unit test" );
    return imm_Y_ref();
}

InfectionOrigin WHInterface::getInfectionOrigin() const{
    if( m_kind == Kind::Vivax || m_kind == Kind::Mock ) return InfectionOrigin::Indigenous;
    return get_infection_origin(infections);
}

// -----  Checkpointing  -----

void WHInterface::checkpoint (istream& stream) {
    numInfs & stream;

    if (numInfs > SlotAllocator::MAX_INFECTIONS)
        throw util::checkpoint_error("numInfs: " + to_string(numInfs));
    
    switch( m_kind ){
    case Kind::Vivax: checkpoint_vivax( stream ); break;
    case Kind::Common: checkpoint_common( stream ); break;
    case Kind::Descriptive: checkpoint_descriptive( stream ); break;
    case Kind::Mock: throw util::unimplemented_exception( "not needed in unit test" );
    }
}
void WHInterface::checkpoint (ostream& stream) {
    numInfs & stream;
    
    switch( m_kind ){
    case Kind::Vivax: checkpoint_vivax( stream ); break;
    case Kind::Common: checkpoint_common( stream ); break;
    case Kind::Descriptive: checkpoint_descriptive( stream ); break;
    case Kind::Mock: throw util::unimplemented_exception( "not needed in unit test" );
    }
}

}
}
