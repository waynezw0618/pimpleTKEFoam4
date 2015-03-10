    // do it parallel
    if (!Pstream::parRun())
    {
        WarningIn(args.executable())
            << "globalIndex class is only useful in parallel code."
            << endl;
    }

    Info<<"searching reference points" <<endl;
    //
    const dictionary& pimpleDict = pimple.dict();
    List<vector> refpts = List<vector>(pimpleDict.lookup("Yrs"));
    labelList ownerProIDlist(refpts.size(),0);
    labelList localIDlist(refpts.size(),0);

    std::vector<volVectorField> RUW(refpts.size(),volVectorField
			     (	
    		        	IOobject
                        (
                            "RUW",
                            runTime.timeName(),
                            mesh,
                            IOobject::READ_IF_PRESENT,
                            IOobject::AUTO_WRITE
                        ),
                		U.component(vector::X)()*fvc::curl(U)
 		        	)
                );
    for(int i=0; i<refpts.size(); i++){
                  RUW[i].rename("RUW"+name(i));
    }
  
    //List< List<vector> > refptsCorr(refpts.size(),List<vector>(Pstream::nProcs(),vector::zero));
    List< List<vector> > refptsCorr(Pstream::nProcs(),List<vector>(refpts.size(),vector::zero));
    List<vector> CorrList(refpts.size(),vector::zero);     
    //Info << "origional refptsCorr are:" << refptsCorr << endl;
    List< List<scalar> > UturbrefListlist(Pstream::nProcs(),List<scalar>(Pstream::nProcs(),0.0)); 
    List<scalar> UturbrefList(refpts.size(),0.0);  

    List< List<scalar> > URMSrefListlist(Pstream::nProcs(),List<scalar>(Pstream::nProcs(),0.0)); 
    List<scalar> URMSrefList(refpts.size(),0.0); 

    forAll(refpts,idx){
        labelList localCellId(Pstream::nProcs(),0);
        localCellId[Pstream::myProcNo()] = mesh.findCell(refpts[idx]);
        reduce(localCellId, maxOp<labelList>());
        
        for (int i=0; i<Pstream::nProcs() ;i++){
            	 if(localCellId[i]!=0){
                       ownerProIDlist[idx]=i;
                       localIDlist[idx]=localCellId[i];
                 }
        }

        if(localIDlist[idx]==0){
                 WarningIn(args.executable())
                 << "the point may not locate in the domain please double check"
                 << endl;
        }
    }
    
    Info <<"referece points are located at Processors :" << ownerProIDlist << nl
         <<"corresponding local ID is : " << localIDlist  <<endl;

    if (Pstream::myProcNo() == 0)
    {
        forAll(ownerProIDlist,idx){
               // we can do this only when first point is not a reference point
               if(ownerProIDlist[idx]==0 && localIDlist[idx]!=0){
                    refptsCorr[0][idx]=mesh.C()[localIDlist[idx]];
               }
        }

        //Info << "after setting at Master, master are:" << refptsCorr << endl;	
        for(label i=1; i<Pstream::nProcs(); i++){
            // create the input stream from processor i
            IPstream vStream(Pstream::blocking, i);
            vStream >> refptsCorr[i];
        }
        // print the list of all vectors on the main proc
        Info << "master is:" << refptsCorr << endl;
    }
    else
    {
        List<vector> localCorrList(refpts.size(),vector::zero);
        forAll(ownerProIDlist,idx){
               if(ownerProIDlist[idx]==Pstream::myProcNo()){
                    localCorrList[idx]=mesh.C()[localIDlist[idx]];
               }
        }
        // print the vector on the processor (so we can compare it later when we print it on the main proc)
        //Sout << "slaver [" << Pstream::myProcNo() << "] localrefCorr = " << localCorrList << endl;

        // create the stream to send to the main proc
        OPstream vectorStream
        (
            Pstream::blocking, 0
        );
        vectorStream << localCorrList;
        
    }

    for(label i=0; i< Pstream::nProcs(); i++){
        for(label j=0; j< refpts.size(); j++){
            //Info << "on Procs [" <<i<<"]th the [" << j 	<<"]th location is at"<< refptsCorr[i][j] << endl;
            if(mag(refptsCorr[i][j])!=0){
                 CorrList[j] = refptsCorr[i][j]; 
            }
        }
    }

    Info << "reference points are located at " << CorrList << endl;


    if (Pstream::myProcNo() == 0)
    {
        forAll(ownerProIDlist,idx){
                // we can do this only when first point is not a reference point
                if(ownerProIDlist[idx]==0 && localIDlist[idx]!=0){
                        URMSrefListlist[0][idx]=URMSMap.component(vector::X)()[localIDlist[idx]];
                }
        }
    
        //Info << "after setting at Master, master are:" << refptsCorr << endl;
        for(label i=1; i<Pstream::nProcs(); i++){
            // create the input stream from processor i
            IPstream vStreamURMS(Pstream::blocking, i);
            vStreamURMS >> URMSrefListlist[i];
        }

    }
    else
    {
        List<scalar> localURMSrefList(refpts.size(),0.0);
        forAll(ownerProIDlist,idx){
                if(ownerProIDlist[idx]==Pstream::myProcNo()){
                      localURMSrefList[idx]=URMSMap.component(vector::X)()[localIDlist[idx]];
                }
        }
    
        OPstream vectorStreamURMS
        (
            Pstream::blocking, 0
        );
        vectorStreamURMS << localURMSrefList;
    }

    for(label i=0; i< Pstream::nProcs(); i++){
        for(label j=0; j< refpts.size(); j++){
            //Info << "on Procs [" <<i<<"]th the [" << j 	<<"]th location is at"<< refptsCorr[i][j] << endl;
            if(mag(refptsCorr[i][j])!=0){
            URMSrefList[j] = URMSrefListlist[i][j];
            }
        }
    }

    Info << "URMS @ master are:" << URMSrefList << endl;
    Pstream::scatter(URMSrefList, Pstream::blocking);
    Pout << "URMS@["<< Pstream::myProcNo() <<"] are :" << URMSrefList << endl;

    forAll(CorrList,idx){
         Info <<"magitude of " << idx <<"th refence point is : " << mag(CorrList[idx]) <<endl;
         if(mag(CorrList[idx])==0 && localIDlist[idx]==0){      //2. mag(CorrList[idx])==0 is wrong!
                FatalErrorIn(args.executable())
                << "refecence point:" << idx << " at " << refpts[idx]
                << "is not at domain please check @_@"
                << abort(FatalError);              
        } 	
    }

