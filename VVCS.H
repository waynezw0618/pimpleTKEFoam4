    // do it parallel
    Info << "calculating VVCS " << endl;
    if (Pstream::myProcNo() == 0)
    {
        forAll(ownerProIDlist,idx){
               if(ownerProIDlist[idx]==0){
                    UturbrefListlist[0][idx]=Uturb.component(vector::X)()[localIDlist[idx]];
               }
        }

        for(label i=1; i<Pstream::nProcs(); i++){
            // create the input stream from processor i
            IPstream vStreamUturb(Pstream::blocking, i);
            vStreamUturb >> UturbrefListlist[i];
        }
 
    }
    else
    {
        List<scalar> localUturbrefList(refpts.size(),0.0);
        forAll(ownerProIDlist,idx){
               if(ownerProIDlist[idx]==Pstream::myProcNo()){
                    localUturbrefList[idx]=Uturb.component(vector::X)()[localIDlist[idx]];
               }
        }        

        // create the stream to send to the main proc
        OPstream vectorStreamUturb
        (
            Pstream::blocking, 0
        );
        vectorStreamUturb << localUturbrefList;

    }


    for(label i=0; i< Pstream::nProcs(); i++){
        for(label j=0; j< refpts.size(); j++){
            if(mag(refptsCorr[i][j])!=0){
                 UturbrefList[j] = UturbrefListlist[i][j];
            }
        }
    } 
    
    Pstream::scatter(UturbrefList, Pstream::blocking);
    Pout<<"Uturref are :" << UturbrefList << endl;

    for(label j=0; j< refpts.size(); j++){
        RUW[j].component(vector::Y) = UturbrefList[j]*Wturb.component(vector::X) \
                                        /(URMSrefList[j]*vorticityRMSMap.component(vector::X));
        RUW[j].component(vector::Y) = UturbrefList[j]*Wturb.component(vector::Y) \
                                        /(URMSrefList[j]*vorticityRMSMap.component(vector::Y));
        RUW[j].component(vector::Y) = UturbrefList[j]*Wturb.component(vector::Z) \
                                        /(URMSrefList[j]*vorticityRMSMap.component(vector::Z));
         if (runTime.outputTime()){
         	 RUW[j].write();
          }
    }      
    //Ruw = dimensionedScalar("one", dimensionSet(0, 1, -1, 0, 0),1.0)*UrefList[0]*fvc::curl(U);
