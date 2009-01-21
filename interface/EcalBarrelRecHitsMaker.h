#ifndef FastSimulation__EcalBarrelRecHitsMaker__h
#define FastSimulation__EcalBarrelRecHitsMaker__h

#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "DataFormats/EcalDigi/interface/EcalDigiCollections.h"
#include "DataFormats/EcalDetId/interface/EcalTrigTowerDetId.h"
#include "FWCore/ParameterSet/interface/InputTag.h"
//#include <boost/cstdint.hpp>

class RandomEngine;
class EcalTrigTowerConstituentsMap;

namespace edm { 
  class ParameterSet;
  class Event;  
  class EventSetup;
}

class EcalBarrelRecHitsMaker
{
 public:
  EcalBarrelRecHitsMaker(edm::ParameterSet const & p, const RandomEngine* );
  ~EcalBarrelRecHitsMaker();

  void loadEcalBarrelRecHits(edm::Event &iEvent, EBRecHitCollection & ecalHits,EBDigiCollection & ecaldigis);
  void init(const edm::EventSetup &es,bool dodigis,bool doMiscalib);

 private:
  void clean();
  void loadPCaloHits(const edm::Event & iEvent);
  void geVtoGainAdc(float e,unsigned& gain,unsigned &adc) const;
  void noisifyTriggerTowers();
  void noisifyTriggerTower(unsigned tthi);
  
 private:
  bool doDigis_;
  bool doMisCalib_;
  double refactor_;
  double refactor_mean_;
  // poor-man Selective Readout
  double threshold_;
  double noise_;
  double calibfactor_;
  const RandomEngine* random_;
  bool noisified_;
  edm::InputTag inputCol_;
  // array (size = 62000) of the energy in the barrel
  std::vector<float> theCalorimeterHits_;
  // array of the hashedindices in the previous array of the cells that received a hit
  std::vector<int> theFiredCells_;
  
  // equivalent of the EcalIntercalibConstants from RecoLocalCalo/EcalRecProducers/src/EcalRecHitProduer.cc
  std::vector<float> theCalibConstants_;

  //array conversion hashedIndex rawId
  std::vector<uint32_t> barrelRawId_;
  float adcToGeV_;
  float geVToAdc1_,geVToAdc2_,geVToAdc3_;
  unsigned minAdc_;
  unsigned maxAdc_;
  float t1_,t2_,sat_;

  const EcalTrigTowerConstituentsMap* eTTmap_;  
  // Array of the DetIds
  std::vector<EcalTrigTowerDetId> theTTDetIds_;
  //Energy of the TT
  std::vector<float> TTEnergy_;
  // shot TTs
  std::vector<unsigned> theFiredTTs_;
  // treated TTs
  std::vector<bool> treatedTTs_;
  // neighboring TT DetIds
  std::vector<std::vector<int> > neighboringTTs_;

  // selective readout threshold
  float SRThreshold_;
  int SREtaSize_;
  int SRPhiSize_;
};

#endif
