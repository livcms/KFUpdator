#include <algorithm>

#include "TrajectorySegmentBuilder.h"

#include "TrackingTools/PatternTools/interface/TempTrajectory.h"

#include "DataFormats/TrajectorySeed/interface/TrajectorySeed.h"
#include "TrackingTools/KalmanUpdators/interface/KFUpdator.h"
#include "TrackingTools/DetLayers/interface/MeasurementEstimator.h"
#include "TrackingTools/PatternTools/interface/Trajectory.h"
#include "TrackingTools/PatternTools/interface/TrajectoryMeasurement.h"
#include "TrackingTools/PatternTools/interface/TrajMeasLessEstim.h"
#include "TrackingTools/MeasurementDet/interface/TrajectoryMeasurementGroup.h"
#include "TrackingTools/DetLayers/interface/DetGroup.h"
#include "TrackingTools/DetLayers/interface/DetLayer.h"
#include "TrajectoryLessByFoundHits.h"
#include "TrackingTools/PatternTools/interface/TrajectoryStateUpdator.h"
#include "TrackingTools/DetLayers/interface/GeomDetCompatibilityChecker.h"
#include "TrackingTools/MeasurementDet/interface/MeasurementDet.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/thread_safety_macros.h"

#include <vector>
#include <typeinfo>
// #define DBG_TSB
// #define STAT_TSB
#include "DataFormats/TrackingRecHit/interface/KfComponentsHolder.h"
#include "DataFormats/Math/interface/ProjectMatrix.h" 
#include <fstream>

#include "TrackingTools/KalmanUpdators/interface/prototype.h"
#include <typeinfo>

namespace {
#ifdef STAT_TSB
  struct StatCount {
    long long totGroup;
    long long totSeg;
    long long totLockHits;
    long long totInvCand;
    long long trunc;
    void zero() { totGroup = totSeg = totLockHits = totInvCand = trunc = 0; }
    void incr(long long g, long long s, long long l) {
      totGroup += g;
      totSeg += s;
      totLockHits += l;
    }
    void truncated() { ++trunc; }
    void invalid() { ++totInvCand; }
    void print() const {
      std::cout << "TrajectorySegmentBuilder stat\nGroup/Seg/Lock/Inv/Trunc " << totGroup << '/' << totSeg << '/'
                << totLockHits << '/' << totInvCand << '/' << trunc << std::endl;
    }
    StatCount() { zero(); }
    ~StatCount() { print(); }
  };
  StatCount statCount;

#else
  struct StatCount {
    void incr(long long, long long, long long) {}
    void truncated() {}
    void invalid() {}
  };
  CMS_THREAD_SAFE StatCount statCount;
#endif

}  // namespace

using namespace std;

TrajectorySegmentBuilder::TempTrajectoryContainer TrajectorySegmentBuilder::segments(const TSOS startingState) {
  //
  // create empty trajectory
  //
  theLockedHits.clear();
  TempTrajectory startingTrajectory(theFullPropagator.propagationDirection(), 0);
  //
  // get measurement groups
  //
  auto&& measGroups =
      theLayerMeasurements->groupedMeasurements(theLayer, startingState, theFullPropagator, theEstimator);

#ifdef DBG_TSB
  cout << "TSB: number of measurement groups = " << measGroups.size() << endl;
  //  theDbgFlg = measGroups.size()>1;
  theDbgFlg = true;
#else
  theDbgFlg = false;
#endif

#ifdef TSB_TRUNCATE
  //  V.I. to me makes things slower...

  //
  // check number of combinations
  //
  constexpr long long MAXCOMB = 100000000;
  long long ncomb(1);
  int ngrp(0);
  bool truncate(false);
  for (auto const& gr : measGroups) {
    ++ngrp;
    int nhit(0);
    for (auto const& m : gr.measurements())
      if
        LIKELY(m.recHitR().isValid()) nhit++;

    if (nhit > 1)
      ncomb *= nhit;
    if
      UNLIKELY(ncomb > MAXCOMB) {
        edm::LogInfo("TrajectorySegmentBuilder")
            << " found " << measGroups.size() << " groups and more than " << static_cast<unsigned int>(MAXCOMB)
            << " combinations - limiting to " << (ngrp - 1) << " groups";
        truncate = true;

        statCount.truncated();

        break;
      }
  }
  //   cout << "Groups / combinations = " << measGroups.size() << " " << ncomb << endl;
  if
    UNLIKELY(truncate && ngrp > 0) measGroups.resize(ngrp - 1);

#endif

#ifdef DBG_TSB
  if (theDbgFlg) {
    int ntot(1);
    for (vector<TMG>::const_iterator ig = measGroups.begin(); ig != measGroups.end(); ++ig) {
      int ngrp(0);
      const vector<TM>& measurements = ig->measurements();
      for (vector<TM>::const_iterator im = measurements.begin(); im != measurements.end(); ++im) {
        if (im->recHit()->isValid())
          ngrp++;
      }
      cout << " " << ngrp;
      if (ngrp > 0)
        ntot *= ngrp;
    }
    cout << endl;
    cout << "TrajectorySegmentBuilder::partialTrajectories:: det ids & hit types / group" << endl;
    for (vector<TMG>::const_iterator ig = measGroups.begin(); ig != measGroups.end(); ++ig) {
      const vector<TM>& measurements = ig->measurements();
      for (vector<TM>::const_iterator im = measurements.begin(); im != measurements.end(); ++im) {
        if (im != measurements.begin())
          cout << " / ";
        if (im->recHit()->det())
          cout << im->recHit()->det()->geographicalId().rawId() << " " << im->recHit()->getType();
        else
          cout << "no det";
      }
      cout << endl;
    }

    //   if ( measGroups.size()>4 ) {
    cout << typeid(theLayer).name() << endl;
    cout << startingState.localError().matrix() << endl;
    //     for (vector<TMG>::const_iterator ig=measGroups.begin();
    // 	 ig!=measGroups.end(); ig++) {
    //       cout << "Nr. of measurements = " << ig->measurements().size() << endl;
    //       const DetGroup& dg = ig->detGroup();
    //       for ( DetGroup::const_iterator id=dg.begin();
    // 	    id!=dg.end(); id++ ) {
    // 	GlobalPoint p(id->det()->position());
    // 	GlobalVector v(id->det()->toGlobal(LocalVector(0.,0.,1.)));
    // 	cout << p.perp() << " " << p.phi() << " " << p.z() << " ; "
    // 	     << v.phi() << " " << v.z() << endl;
    //       }
    //     }
    //   }
  }
#endif

  TempTrajectoryContainer candidates = addGroup(startingTrajectory, measGroups.begin(), measGroups.end());

  if
    UNLIKELY(theDbgFlg) cout << "TSB: back with " << candidates.size() << " candidates" << endl;

  //
  // add invalid hit - try to get first detector hit by the extrapolation
  //

  updateWithInvalidHit(startingTrajectory, measGroups, candidates);

  if
    UNLIKELY(theDbgFlg) cout << "TSB: " << candidates.size() << " candidates after invalid hit" << endl;

  statCount.incr(measGroups.size(), candidates.size(), theLockedHits.size());

  theLockedHits.clear();

  return candidates;
}



template <typename T>
void writeVectorToFile (T t, ofstream& file){
  for (auto j = t.begin(); j != t.end(); j++){
            file.write(reinterpret_cast<char*> (& *j), sizeof(*j));
  //return 0;
        }
}


void TrajectorySegmentBuilder::updateTrajectory(TempTrajectory& traj, TM tm) const {
  auto&& predictedState = tm.predictedState();
  auto&& hit = tm.recHit();
//  
  if( hit -> isValid()){
    auto && x =  predictedState.localParameters().vector();
    auto && C = predictedState.localError().matrix();

    ProjectMatrix<double, 5, 2> pf;

    AlgebraicROOTObject<2>::Vector r, rMeas; 
 
    using ROOT::Math::SMatrixNoInit;
    AlgebraicROOTObject<2, 2>::SymMatrix V(SMatrixNoInit{}), VMeas(SMatrixNoInit{});

    KfComponentsHolder holder;
    holder.template setup<2>(&r, &V, &pf, &rMeas, &VMeas, x, C);
    const TrackingRecHit& aRecHit = *hit; 
    aRecHit.getKfComponents(holder);
  
    std::vector<double> xvec(x.begin(), x.end());
    std::vector<double> Cvec(C.begin(), C.end()); 
    std::vector<double> rvec(r.begin(), r.end()); 
    std::vector<double> Vvec(V.begin(), V.end()); 
    std::vector<double> VMeasvec(VMeas.begin(), VMeas.end());   

   static int printonce= 0; 
    if(printonce < 5){
      std::cout<<" x is "<<x<<std::endl; 
      std::cout<<" C is "<<C<<std::endl; 
      std::cout<<" r is "<<r<<std::endl; 
      std::cout<<" rMeas is "<<rMeas<<std::endl; 
      std::cout<<" V is "<<V<<std::endl; 
      std::cout<<" VMeas is "<<VMeas<<std::endl; 
      for (int i = 0; i < 5; i++){ 
		std::cout<<"x written to file "<< xvec[i]<<std::endl; 
}
      printonce++; 

}
    //if (printonce < 5){
//	for (int i = 0; i < 5; i++){ 
//		std::cout<<"x written to file "<< xvec[i]<<std::endl; 
   // }
   //}


     ofstream wf("KFUpdatorInputsPU200_numHit10.dat", ios::binary| std::ios::ate |std::ios::app); 
     if(wf){
      writeVectorToFile(xvec, wf); 
     // double numMeasurements = measurements.size(); 
      writeVectorToFile(Cvec, wf); 
      writeVectorToFile(rvec, wf); 
      //writeVectorToFile(rMeas, wf);
      writeVectorToFile(Vvec, wf);  
      writeVectorToFile(VMeasvec, wf);   
     wf.close(); 
}
    output temp = testfunction(xvec, Cvec, rvec, Vvec, VMeasvec); 
    std::vector tempfsv = temp.finalStateVector; 
    std::vector tempfse  = temp.finalStateError; 

    static bool checkprint = true; 
    if (checkprint){
    for( auto i = tempfsv.begin(); i<tempfsv.end(); i++){ 
	std::cout<<"tempfsv is "<<*i<<std::endl; 
    }
    for( auto i = tempfse.begin(); i<tempfse.end(); i++){ 
	std::cout<<"tempfse is "<<*i<<std::endl; 
    }
    checkprint = false; 
    }

ofstream wf1("KFUpdatorCorrectOutputPU200_numHits10.dat", ios::binary| std::ios::ate |std::ios::app); 
    if(wf1){
      writeVectorToFile(tempfsv, wf1); 
     // double numMeasurements = measurements.size(); 
      writeVectorToFile(tempfse, wf1); 
     wf1.close(); 
	}
 // AlgebraicVector5 newtestfsv = {tempfsv[0], tempfsv[1], tempfsv[2], tempfsv[3], tempfsv[4]}; 
    //AlgebraicSymMatrix55 newtestfs = {tempfse[0], tempfse[5], tempfse[6], tempfse[10], tempfse[11], tempfse[12], tempfse[15], tempfse[16], tempfse[17], tempfse[18], tempfse[20], tempfse[21], tempfse[22], tempfse[23], tempfse[24]};
  // ROOT::Math::SVector<double, 15> v(); 
   //AlgebraicSymMatrix55 newtestfse(tempfse.begin(), tempfse.end());
//   auto&& upState  = TrajectoryStateOnSurface(LocalTrajectoryParameters(newtestfsv, predictedState.localParameters().pzSign()),
//                                      LocalTrajectoryError(newtestfse),
//                                      predictedState.surface(),
//                                      &(predictedState.globalParameters().magneticField()),
//                                      predictedState.surfaceSide());
//  
//    if (hit->isValid()) {
    auto&& upState = theUpdator.update(predictedState, *hit);
    static int printupstate=0; 
    if(printupstate  < 5 ){
    std::cout<<"the upstate is "<<upState<<std::endl; 
    //  std::cout<<"local trajectory error test "<<LocalTrajectoryError()<<std::endl;  
      //std::cout<<" type of upstate is "<<typeid(upState).name()<<std::endl;
      printupstate++; 
    }
    traj.emplace(std::move(predictedState), std::move(upState), std::move(hit), tm.estimate(), tm.layer());
 
  //       TrajectoryMeasurement tm(traj.lastMeasurement());
    //    if ( tm.updatedState().isValid() ) {
      //     if ( !hit.det().surface()->bounds().inside(tm.updatedState().localPosition(),
   // static bool printhit = true; 
    //if(printhit){

    //std::cout<<"hit at "<<hit.localPosition()<< " : " << endl;
   // printonce = false; 
   // }
   
   //  	tm.updatedState().localError().positionError(),3.f) ) {

     //	cout << "Incompatibility after update for det at " << hit.det().position() << ":" << endl;
     //	cout << tm.predictedState().localPosition() << " "
     //	     << tm.predictedState().localError().positionError() << endl;
     //	cout << hit.localPosition() << " " << hit.localPositionError() << endl;
     //	cout << tm.updatedState().localPosition() << " "
     //	     << tm.updatedState().localError().positionError() << endl;
    //       }
    //     }
  } else {
    traj.emplace(std::move(predictedState), std::move(hit), 0, tm.layer());
  }

}

TrajectorySegmentBuilder::TempTrajectoryContainer TrajectorySegmentBuilder::addGroup(TempTrajectory const& traj,
                                                                                     vector<TMG>::const_iterator begin,
                                                                                     vector<TMG>::const_iterator end) {
  vector<TempTrajectory> ret;
  if (begin == end) {
    //std::cout << "TrajectorySegmentBuilder::addGroup" << " traj.empty()=" << traj.empty() << "EMPTY" << std::endl;
    if
      UNLIKELY(theDbgFlg) cout << "TSB::addGroup : no groups left" << endl;
    if (!traj.empty())
      ret.push_back(traj);
    return ret;
  }

  if
    UNLIKELY(theDbgFlg)
  cout << "TSB::addGroup : traj.size() = " << traj.measurements().size() << " first group at "
       << &(*begin)
       //        << " nr. of candidates = " << candidates.size()
       << endl;
  

  TempTrajectoryContainer updatedTrajectories;
  updatedTrajectories.reserve(2);
  if (traj.measurements().empty()) {
    if (theMaxCand == 1) {
      auto&& firstMeasurements = unlockedMeasurements(begin->measurements());
      if (!firstMeasurements.empty())
        updateCandidatesWithBestHit(traj, std::move(firstMeasurements.front()), updatedTrajectories);
    } else {
      updateCandidates(traj, begin->measurements(), updatedTrajectories);
    }
    if
      UNLIKELY(theDbgFlg)
    cout << "TSB::addGroup : updating with first group - " << updatedTrajectories.size() << " trajectories" << endl;
  } else {
    auto&& meas = redoMeasurements(traj, begin->detGroup());
    if (!meas.empty()) {
      if (theBestHitOnly) {
        updateCandidatesWithBestHit(traj, std::move(meas.front()), updatedTrajectories);
      } else {
        updateCandidates(traj, std::move(meas), updatedTrajectories);
      }
      if
        UNLIKELY(theDbgFlg)
      cout << "TSB::addGroup : updating" << updatedTrajectories.size() << " trajectories-1" << endl;
    }
  }
  // keep old trajectory
  //
  updatedTrajectories.push_back(traj);

  if (begin + 1 != end) {
    ret.reserve(4);  // a good upper bound
    for (auto const& ut : updatedTrajectories) {
      if
        UNLIKELY(theDbgFlg)
      cout << "TSB::addGroup : trying to extend candidate at " << &ut << " size " << ut.measurements().size() << endl;
      vector<TempTrajectory>&& finalTrajectories = addGroup(ut, begin + 1, end);
      if
        UNLIKELY(theDbgFlg)
      cout << "TSB::addGroup : " << finalTrajectories.size() << " finalised candidates before cleaning" << endl;
      //B.M. to be ported later
      // V.I. only mark invalidate
      cleanCandidates(finalTrajectories);

      if
        UNLIKELY(theDbgFlg) {
          int ntf = 0;
          for (auto const& t : finalTrajectories)
            if (t.isValid())
              ++ntf;
          cout << "TSB::addGroup : got " << ntf << " finalised candidates" << endl;
        }

      for (auto& t : finalTrajectories)
        if (t.isValid())
          ret.push_back(std::move(t));

      //        ret.insert(ret.end(),make_move_iterator(finalTrajectories.begin()),
      //	   make_move_iterator(finalTrajectories.end()));
    }
  } else {
    ret.reserve(updatedTrajectories.size());
    for (auto& t : updatedTrajectories)
      if (!t.empty())
        ret.push_back(std::move(t));
  }

  //std::cout << "TrajectorySegmentBuilder::addGroup" <<
  //             " traj.empty()=" << traj.empty() <<
  //             " end-begin=" << (end-begin)  <<
  //             " #updated=" << updatedTrajectories.size() <<
  //             " #result=" << ret.size() << std::endl;
  return ret;
}

//ofstream wf("testprint.dat", ios::binary| std::ios::ate |std::ios::app); 
//ofstream wf("testprint.dat", ios::binary, fstream::app); 
//fstream wf; 
//wf.open("testprint.dat", fstream::app)

void TrajectorySegmentBuilder::updateCandidates(TempTrajectory const& traj,
                                                const vector<TM>& measurements,
                                                TempTrajectoryContainer& candidates) {
  //
  // generate updated candidates with all valid hits
  
  for (auto im = measurements.begin(); im != measurements.end(); ++im) {
    if (im->recHit()->isValid()) {
      
      candidates.push_back(traj); 
	
// Liv's code for writing to file! 
//	TM tm = *im; 
//        auto && predictedState = tm.predictedState(); 
//        auto && hit = tm.recHit();
////	if( hit -> isValid()){
//
//	  auto && x =  predictedState.localParameters().vector();
//  	  auto && C = predictedState.localError().matrix();
//  	  ProjectMatrix<double, 5, 2> pf;
//
//  	  AlgebraicROOTObject<2>::Vector r, rMeas;
//
//  	  using ROOT::Math::SMatrixNoInit;
//  	AlgebraicROOTObject<2, 2>::SymMatrix V(SMatrixNoInit{}), VMeas(SMatrixNoInit{});
//
//  	KfComponentsHolder holder;
//  	holder.template setup<2>(&r, &V, &pf, &rMeas, &VMeas, x, C);
//  	const TrackingRecHit& aRecHit = *hit;
//  	aRecHit.getKfComponents(holder);
////        printcand++; 
//
//       ofstream wf("KFUpdatorInputsPU0_numHit10.dat", ios::binary| std::ios::ate |std::ios::app); 
//       if  (wf){
//         writeVectorToFile(x, wf); 
//        // double numMeasurements = measurements.size(); 
//         writeVectorToFile(C, wf); 
//         writeVectorToFile(r, wf); 
//	 //writeVectorToFile(rMeas, wf);
//         writeVectorToFile(V, wf);  
//	 writeVectorToFile(VMeas, wf);   
//        wf.close(); 
//}
// End of Liv's code for writing to file   
      
      updateTrajectory(candidates.back(), *im);
      if (theLockHits)
        lockMeasurement(*im);

}
    }
//  wf.cloie();
}



void TrajectorySegmentBuilder::updateCandidatesWithBestHit(TempTrajectory const& traj,
                                                           TM measurement,
                                                           TempTrajectoryContainer& candidates) {
  // here we arrive with only valid hits and sorted.
  //so the best is the first!

  if (theLockHits)
    lockMeasurement(measurement);
  candidates.push_back(traj); 
  
//	TM tm = measurement; 
 //       auto && predictedState = tm.predictedState(); 
//        auto && hit = tm.recHit();
//	if( hit -> isValid()){
//
//	  auto && x =  predictedState.localParameters().vector();
//  	  auto && C = predictedState.localError().matrix();
//  	  ProjectMatrix<double, 5, 2> pf;
//
//  	  AlgebraicROOTObject<2>::Vector r, rMeas;
//
//  	  using ROOT::Math::SMatrixNoInit;
//  	AlgebraicROOTObject<2, 2>::SymMatrix V(SMatrixNoInit{}), VMeas(SMatrixNoInit{});
//
//  	KfComponentsHolder holder;
//  	holder.template setup<2>(&r, &V, &pf, &rMeas, &VMeas, x, C);
//  	const TrackingRecHit& aRecHit = *hit;
//  	aRecHit.getKfComponents(holder);
////        printcand++; 
// 
//       ofstream wf("testprint.dat", ios::binary| std::ios::ate |std::ios::app); 
//       if  (wf){
//         writeVectorToFile(x, wf); 
//        // double numMeasurements = measurements.size(); 
//         writeVectorToFile(C, wf); 
//         //wf.write(reinterpret_cast<char*> (& numMeasurements), sizeof(numMeasurements)); 
//         //wf.write( reinterpret_cast<char*>(& rprint), sizeof(rprint)); 
//         writeVectorToFile(r, wf); 
//	 writeVectorToFile(rMeas, wf);
//         writeVectorToFile(V, wf);  
//	 writeVectorToFile(VMeas, wf);   
//	//wf.write(reinterpret_cast<char*> (& C), sizeof(C)); 
//	//      
//	//        struct printvec{ 
////          //int nhits 
////          double rvec1, rvec2; 
////       } rprint; 
////        
////        rprint.rvec1  = r[0];
////        rprint.rvec2 = r[1]; 
//        wf.close(); 
        //file.open("testprint.dat", ios::trunc); 
        //file<< r[0] << r[1]; 
        //file.close(); 
//	       }                
//}
   
//	std::cout << "number of hits found is " << measurements.size() <<endl; 
	//printcand = false; 
  //    }
  updateTrajectory(candidates.back(), std::move(measurement));
}

vector<TrajectoryMeasurement> TrajectorySegmentBuilder::redoMeasurements(const TempTrajectory& traj,
                                                                         const DetGroup& detGroup) const {
  vector<TM> result;
  //
  // loop over all dets
  //
  if
    UNLIKELY(theDbgFlg) cout << "TSB::redoMeasurements : nr. of measurements / group =";

  tracking::TempMeasurements tmps;

  for (auto const& det : detGroup) {
    pair<bool, TrajectoryStateOnSurface> compat = GeomDetCompatibilityChecker().isCompatible(
        det.det(), traj.lastMeasurement().updatedState(), theGeomPropagator, theEstimator);

    if
      UNLIKELY(theDbgFlg && !compat.first) std::cout << " 0";

    if (!compat.first)
      continue;

    MeasurementDetWithData mdet = theLayerMeasurements->idToDet(det.det()->geographicalId());
    // verify also that first (and only!) not be inactive..
    if (mdet.measurements(compat.second, theEstimator, tmps) && tmps.hits[0]->isValid())
      for (std::size_t i = 0; i != tmps.size(); ++i)
        result.emplace_back(compat.second, std::move(tmps.hits[i]), tmps.distances[i], &theLayer);

    if
      UNLIKELY(theDbgFlg) std::cout << " " << tmps.size();
    tmps.clear();
  }

  if
    UNLIKELY(theDbgFlg) cout << endl;

  std::sort(result.begin(), result.end(), TrajMeasLessEstim());

  return result;
}

void TrajectorySegmentBuilder::updateWithInvalidHit(TempTrajectory& traj,
                                                    const vector<TMG>& groups,
                                                    TempTrajectoryContainer& candidates) const {
  //
  // first try to find an inactive hit with dets crossed by the prediction,
  //   then take any inactive hit
  //
  // loop over groups
  for (int iteration = 0; iteration < 2; iteration++) {
    for (auto const& gr : groups) {
      auto const& measurements = gr.measurements();
      for (auto im = measurements.rbegin(); im != measurements.rend(); ++im) {
        auto const& hit = im->recHitR();
        if ((hit.getType() == TrackingRecHit::valid) | (hit.getType() == TrackingRecHit::missing))
          continue;
        //
        // check, if the extrapolation traverses the Det or
        // if 2nd iteration
        //
        if (hit.det()) {
          auto const& predState = im->predictedState();
          if (iteration > 0 ||
              (predState.isValid() && hit.det()->surface().bounds().inside(predState.localPosition()))) {
            // add the hit
            candidates.push_back(traj);
            updateTrajectory(candidates.back(), *im);
            if
              UNLIKELY(theDbgFlg)
            cout << "TrajectorySegmentBuilder::updateWithInvalidHit "
                 << "added inactive hit" << endl;
            return;
          }
        }
      }
    }
  }
  //
  // No suitable inactive hit: add a missing one
  //
  for (int iteration = 0; iteration < 2; iteration++) {
    //
    // loop over groups
    //
    for (auto const& gr : groups) {
      auto const& measurements = gr.measurements();
      for (auto im = measurements.rbegin(); im != measurements.rend(); ++im) {
        // only use invalid hits
        auto const& hit = im->recHitR();
        if
          LIKELY(hit.isValid()) continue;

        // check, if the extrapolation traverses the Det
        auto const& predState = im->predictedState();
        if (iteration > 0 || (predState.isValid() && hit.surface()->bounds().inside(predState.localPosition()))) {
          // add invalid hit
          candidates.push_back(traj);
          updateTrajectory(candidates.back(), *im);
          return;
        }
      }
    }
    if
      UNLIKELY(theDbgFlg && iteration == 0)
    cout << "TrajectorySegmentBuilder::updateWithInvalidHit: "
         << " did not find invalid hit on 1st iteration" << endl;
  }

  if
    UNLIKELY(theDbgFlg)
  cout << "TrajectorySegmentBuilder::updateWithInvalidHit: "
       << " did not find invalid hit" << endl;
}

vector<TrajectoryMeasurement> TrajectorySegmentBuilder::unlockedMeasurements(const vector<TM>& measurements) const {
  //   if ( !theLockHits )  return measurements;

  vector<TM> result;
  result.reserve(measurements.size());

  //RecHitEqualByChannels recHitEqual(false,true);

  for (auto const& m : measurements) {
    auto const& testHit = m.recHitR();
    if
      UNLIKELY(!testHit.isValid()) continue;
    bool found(false);
    if
      LIKELY(theLockHits) {
        for (auto const& h : theLockedHits) {
          if (h->hit()->sharesInput(testHit.hit(), TrackingRecHit::all)) {
            found = true;
            break;
          }
        }
      }
    if
      LIKELY(!found) result.push_back(m);
  }
  return result;
}

void TrajectorySegmentBuilder::lockMeasurement(const TM& measurement) { theLockedHits.push_back(measurement.recHit()); }

// ================= B.M. to be ported later ===============================
void TrajectorySegmentBuilder::cleanCandidates(vector<TempTrajectory>& candidates) const {
  //
  // remove candidates which are subsets of others
  // assumptions: no invalid hits and no duplicates
  //
  if (candidates.size() <= 1)
    return;
  //RecHitEqualByChannels recHitEqual(false,true);
  //
  const int NC = candidates.size();
  int index[NC];
  for (int i = 0; i != NC; ++i)
    index[i] = i;
  std::sort(index, index + NC, [&candidates](int i, int j) { return lessByFoundHits(candidates[i], candidates[j]); });
  //   cout << "SortedCandidates.foundHits";
  //   for (auto i1 : index)
  //     cout << " " << candidates[i1].foundHits();
  //   cout << endl;
  //
  for (auto i1 = index; i1 != index + NC - 1; ++i1) {
    // get measurements of candidate to be checked
    const TempTrajectory::DataContainer& measurements1 = candidates[*i1].measurements();
    for (auto i2 = i1 + 1; i2 != index + NC; ++i2) {
      // no duplicates: two candidates of same size are different
      if (candidates[*i2].foundHits() == candidates[*i1].foundHits())
        continue;
      // get measurements of "reference"
      const TempTrajectory::DataContainer& measurements2 = candidates[*i2].measurements();
      //
      // use the fact that TMs are ordered:
      // start search in trajectory#1 from last hit match found
      //
      bool allFound(true);
      TempTrajectory::DataContainer::const_iterator from2 = measurements2.rbegin(), im2end = measurements2.rend();
      for (TempTrajectory::DataContainer::const_iterator im1 = measurements1.rbegin(), im1end = measurements1.rend();
           im1 != im1end;
           --im1) {
        // redundant protection - segments should not contain invalid RecHits
        // assert( im1->recHit()->isValid());
        bool found(false);
        for (TempTrajectory::DataContainer::const_iterator im2 = from2; im2 != im2end; --im2) {
          // redundant protection - segments should not contain invalid RecHits
          // assert (im2->recHit()->isValid());
          if (im1->recHitR().hit()->sharesInput(im2->recHitR().hit(), TrackingRecHit::all)) {
            found = true;
            from2 = im2;
            --from2;
            break;
          }
        }
        if (!found) {
          allFound = false;
          break;
        }
      }
      if (allFound) {
        candidates[*i1].invalidate();
        statCount.invalid();
      }
    }
  }
}

//==================================================
