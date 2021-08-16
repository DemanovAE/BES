//***********************CALCULATE Q-VECTOR**********************//
TVector2 CalculateQvector(Int_t harmonic, StFemtoTrack *const &track){
  
  TVector2 qv(0.,0.);
  qv.Set( TMath::Cos( harmonic * track -> phi() ) , TMath::Sin( harmonic * track -> phi() ) );
  return qv;

}// CalculateQvector(){}

Int_t GetEtaDirection(StFemtoTrack *const &track){

  if(track -> eta() < 0.) return 0;
  if(track -> eta() > 0.) return 1;
  return -1;
}

Int_t GetBinVtxZ(StFemtoEvent *const &event, const Int_t _energy){
  
  TVector3 pVtx = event->primaryVertex();

  if(pVtx.Z() == (-1.0 * CutVtxZ.at(_energy))) return 0;
  if(pVtx.Z() == CutVtxZ.at(_energy)) return (nBinVtxZ_PID-1);

  Int_t bin = -1;

  bin = (Int_t)( TMath::Abs(CutVtxZ.at(_energy) + pVtx.Z()) / (2.0*CutVtxZ.at(_energy) / nBinVtxZ_PID));
  
  //std::cout << bin <<"\t\t" << pVtx.Z()<<std::endl;

  if(bin > nBinVtxZ_PID)return -1;

  return bin;

}

Int_t GetBinEta(StFemtoTrack *const &track){

  for(Int_t i = nEtaGap-1 ; i >= 0; i--){
    if( TMath::Abs(track -> eta()) > EtaVecPID[i]){ 
      return i;
    } 
  }
  return -1;
}

Float_t GetMass(Int_t particle){
  if(particle==0){
    return pion_mass;
  }
  if(particle==1){
    return kaon_mass;
  }
  if(particle==2){
    return proton_mass;
  }
  return 100.;
}

Float_t GetRapidity(StFemtoTrack *const &track, Int_t particle){

  Float_t E = sqrt( pow( track->p(), 2) + pow( GetMass(particle) ,2) );
  return 0.5 * log( (E + track->pMom().Z() ) / ( E - track->pMom().Z() ) );
}

Double_t GetWeight(StFemtoTrack *const &track){
  
  Double_t w;
  if (track->pt() < 2.0){
    w = track->pt();
  }
  else{
    w = 2.0;
  }
  return w;
}

int GetCharge(StFemtoTrack *const &track){

  if((Int_t)track -> charge() > 0) return 0;
  if((Int_t)track -> charge() < 0) return 1;

  return -1;
}

Int_t GetBinPtRange(StFemtoTrack *const &track){

  if ( track->isTofTrack() ){
    for(Int_t i = 0; i < (int)ptBinRange.size()-1; i++){
      if( track->pt() >= ptBinRange[i] && track->pt() < ptBinRange[i+1]){
        return i;
      }
    }
  }
  return -1;
}

Bool_t TofMatchedCut(StFemtoDst *const &dst, Int_t cutTofMatched){

    Int_t nTrack = dst->numberOfTracks();
    Int_t number_tof=0;

    for(Int_t iTrk=0; iTrk<nTrack; iTrk++) {
      StFemtoTrack *femtoTrack = dst->track(iTrk);
      if ( !femtoTrack ) continue;
      if ( femtoTrack->isTofTrack()){
        number_tof++;
      }
      if(number_tof > cutTofMatched) return true;
    }

    return false;

}