/* This file is part of OpenMalaria.
 * 
 * Copyright (C) 2005-2026 Swiss Tropical and Public Health Institute
 * Copyright (C) 2005-2015 Liverpool School Of Tropical Medicine
 * Copyright (C) 2020-2026 University of Basel
 * Copyright (C) 2025-2026 The Kids Research Institute Australia
 *
 * OpenMalaria is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#include <iostream>
#include <fstream>
#include <gzstream/gzstream.h>

#include "Global.h"
#include "Model.h"
#include "mon/Continuous.h"
#include "mon/Monitoring.h"
#include "interventions/Interventions.h"
#include "Host/NeonatalMortality.h"
#include "Clinical/ClinicalModel.h"

#include "checkpoint.h"
#include "Transmission/TransmissionModel.h"

#include "schema/scenario.h"
#include "Transmission/VectorModel.h"
#include "Transmission/NonVectorModel.h"
#include "Transmission/Anopheles/SimpleMPDAnophelesModel.h"

// This field order is the existing checkpoint format, not the in-memory layout.
namespace OM::Transmission {
template<class S>
void TransmissionModel::checkpointFields(S& stream) {
    simulationMode &stream;
    interventionMode &stream;
    initialisationEIR &stream;
    laggedKappa &stream;
    annualEIR &stream;
    _annualAverageKappa &stream;
    _sumAnnualKappa &stream;
    tsAdultEIR &stream;
    surveyInputEIR &stream;
    surveySimulatedEIR &stream;
    surveySimulatedEIR_i &stream;
    surveySimulatedEIR_l &stream;
    lastSurveyTime &stream;
    adultAge &stream;
    numTransmittingHumans &stream;
}
void TransmissionModel::checkpoint(istream& stream) { checkpointFields(stream); }
void TransmissionModel::checkpoint(ostream& stream) { checkpointFields(stream); }

template<class S>
void NonVectorModel::checkpointFields(S& stream) {
    TransmissionModel::checkpoint(stream);
    nSpore &stream;
    interventionEIR &stream;
    initialKappa &stream;
}
void NonVectorModel::checkpoint(istream& stream) { checkpointFields(stream); }
void NonVectorModel::checkpoint(ostream& stream) { checkpointFields(stream); }

template<class S>
void VectorModel::checkpointFields(S& stream) {
    TransmissionModel::checkpoint(stream);
    m_rng.checkpoint(stream);
    initIterations &stream;
    for (auto &s : species)
        s->checkpoint(stream);
}
void VectorModel::checkpoint(istream& stream) { checkpointFields(stream); }
void VectorModel::checkpoint(ostream& stream) { checkpointFields(stream); }

template<class S>
void Anopheles::AnophelesModel::checkpointFields(S& stream) {
    mosq.seekingDeathRate & stream;
    mosq.seekingDuration & stream;
    mosq.probMosqSurvivalOvipositing & stream;
    seekingDeathRateIntervs & stream;
    probDeathOvipositingIntervs & stream;
    partialEIR_i & stream;
    partialEIR_l & stream;
    EIRRotateAngle & stream;
    FSCoeffic & stream;
    forcedS_v & stream;
    initNvFromSv & stream;
    initOvFromSv & stream;
    emergenceReduction & stream;
    mosqEmergeRate & stream;
    quinquennialS_v & stream;
    initNv0FromSv & stream;
    initSvFromEIR & stream;
    mosq.restDuration & stream;
    mosq.EIPDuration & stream;
    N_v_length & stream;
    P_A & stream;
    P_df & stream;
    P_dif_i & stream;
    P_dif_l & stream;
    P_dff & stream;
    N_v & stream;
    O_v_i & stream;
    O_v_l & stream;
    S_v_i & stream;
    S_v_l & stream;
    P_Amu & stream;
    P_A1 & stream;
    P_Ah & stream;
    fArray & stream;
    ftauArray & stream;
    uninfected_v & stream;
    timeStep_N_v0 & stream;

}
void Anopheles::AnophelesModel::checkpoint(istream& stream) { checkpointFields(stream); }
void Anopheles::AnophelesModel::checkpoint(ostream& stream) { checkpointFields(stream); }

template<class S>
void Anopheles::SimpleMPDAnophelesModel::checkpointFields(S& stream) {
    AnophelesModel::checkpoint(stream);

    quinquennialOvipositing &stream;
    developmentDuration &stream;
    probPreadultSurvival &stream;
    fEggsLaidByOviposit &stream;
    invLarvalResources &stream;
    nOvipositingDelayed &stream;

}
void Anopheles::SimpleMPDAnophelesModel::checkpoint(istream& stream) { checkpointFields(stream); }
void Anopheles::SimpleMPDAnophelesModel::checkpoint(ostream& stream) { checkpointFields(stream); }
}

namespace OM
{

/** @brief checkpointing functions
*
* readCheckpoint/writeCheckpoint prepare to read/write the file,
* and read/write read and write the actual data. */
//@{
int readCheckpointNum (const string &checkpointFileName)
{
    ifstream checkpointFile;
    checkpointFile.open(checkpointFileName, fstream::in);
    int checkpointNum=0;
    checkpointFile >> checkpointNum;
    checkpointFile.close();
    if (!checkpointFile)
        throw util::checkpoint_error ("error reading from file \"checkpoint\"");
    return checkpointNum;
}

void checkpoint (istream& stream, SimTime &endTime, SimTime &estEndTime, Model& model)
{
    try
    {
        util::checkpoint::header (stream);
        util::CommandLine::staticCheckpoint (stream);
        model.population.neonatalMortality & stream;
        model.reporting.newInfections & stream;
        model.reporting.infantMortality.deaths & stream;
        model.reporting.infantMortality.intervalsAtRisk & stream;
        WithinHost::Genotypes::staticCheckpoint(stream);
        mon::Continuous::checkpoint(model.reporting.continuous, stream);
        mon::checkpoint( stream );
#       ifdef OM_STREAM_VALIDATOR
        util::StreamValidator & stream;
#       endif
        
        sim::s_interv & stream;
        endTime & stream;
        estEndTime & stream;
        *model.transmission & stream;
        checkpoint(model.population, stream, model.parameters, model.masterRng, model.interventions.humanComponents);
        model.interventions.checkpoint(stream);
        model.interventions.loadFromCheckpoint(model.population.humans, *model.transmission);
        
        // read last, because other loads may use random numbers or expect time
        // to be negative
        sim::s_t0 & stream;
        sim::s_t1 & stream;
        model.masterRng.checkpoint(stream);
    } catch (const util::checkpoint_error& e) { // append " (pos X of Y bytes)"
        ostringstream pos;
        pos<<" (pos "<<stream.tellg()<<" of ";
        stream.ignore (numeric_limits<streamsize>::max()-1);    // skip to end of file
        pos<<stream.tellg()<<" bytes)";
        throw util::checkpoint_error( e.what() + pos.str() );
    }
    
    
    stream.ignore (numeric_limits<streamsize>::max()-1);        // skip to end of file
    if (stream.gcount () != 0) {
        ostringstream msg;
        msg << "Checkpointing file has " << stream.gcount() << " bytes remaining." << endl;
        throw util::checkpoint_error (msg.str());
    } else if (stream.fail())
        throw util::checkpoint_error ("stream read error");
}

void checkpoint (ostream& stream, SimTime &endTime, SimTime &estEndTime, Model& model) {
    util::checkpoint::header (stream);
    if (!stream.good())
        throw util::checkpoint_error ("Unable to write to file");

    util::CommandLine::staticCheckpoint (stream);
    model.population.neonatalMortality & stream;
    model.reporting.newInfections & stream;
    model.reporting.infantMortality.deaths & stream;
    model.reporting.infantMortality.intervalsAtRisk & stream;
    WithinHost::Genotypes::staticCheckpoint(stream);
    mon::Continuous::checkpoint(model.reporting.continuous, stream);
    mon::checkpoint( stream );
# ifdef OM_STREAM_VALIDATOR
    util::StreamValidator & stream;
# endif
    
    sim::s_interv & stream;
    endTime & stream;
    estEndTime & stream;
    *model.transmission & stream;
    checkpoint(model.population, stream);
    model.interventions.checkpoint( stream );
    
    sim::s_t0 & stream;
    sim::s_t1 & stream;
    model.masterRng.checkpoint(stream);
    
    if (stream.fail())
        throw util::checkpoint_error ("stream write error");
}

void writeCheckpoint(const bool startedFromCheckpoint, const string &checkpointFileName, SimTime &endTime, SimTime &estEndTime, Model& model)
{
    // We alternate between two checkpoints, in case program is closed while writing.
    const int NUM_CHECKPOINTS = 2;
    
    int oldCheckpointNum = 0, checkpointNum = 0;
    if (startedFromCheckpoint)
    {
        oldCheckpointNum = readCheckpointNum(checkpointFileName);
        checkpointNum = mod_nn(oldCheckpointNum + 1, NUM_CHECKPOINTS); // Get next checkpoint number:
    }
    
    {   // Open the next checkpoint file for writing:
        ostringstream name;
        name << checkpointFileName << checkpointNum << ".gz";
        ogzstream out(name.str().c_str(), ios::out | ios::binary);
        checkpoint (out, endTime, estEndTime, model);
        out.close();
    }
    
    {   // Indicate which is the latest checkpoint file.
        ofstream checkpointFile;
        checkpointFile.open(checkpointFileName,ios::out);
        checkpointFile << checkpointNum;
        checkpointFile.close();
        if (!checkpointFile)
            throw util::checkpoint_error ("error writing to file \"checkpoint\"");
    }

    // Truncate the old checkpoint to save disk space, when it existed
    if( oldCheckpointNum != checkpointNum ){
        ostringstream name;
        name << checkpointFileName << oldCheckpointNum << ".gz";
        ofstream out(name.str().c_str(), ios::out | ios::binary);
        out.close();
    }
}

void readCheckpoint(const string &checkpointFileName, SimTime &endTime, SimTime &estEndTime, Model& model)
{
    int checkpointNum = readCheckpointNum(checkpointFileName);
    
    // Open the latest file
    ostringstream name;
    name << checkpointFileName << checkpointNum << ".gz";
    igzstream in(name.str().c_str(), ios::in | ios::binary);
    //Note: gzstreams are considered "good" when file not open!
    if ( !( in.good() && in.rdbuf()->is_open() ) )
        throw util::checkpoint_error ("Unable to read file");
    checkpoint (in, endTime, estEndTime, model);
    in.close();
  
    cerr << sim::inSteps(sim::now()) << "t loaded checkpoint" << endl;
}

}
