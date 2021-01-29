// make sure the premixed flame front always in the center of the domain, only for regular hexahedron mesh
// only for serial simulation and one-dimensional case now
{
	Info<< "Reading moveFields properties\n" << endl;

	IOdictionary moveFieldsProperties
	(
		IOobject
		(
			"moveFieldsProperties",
			runTime.constant(),
			mesh,
			IOobject::MUST_READ,
			IOobject::NO_WRITE,
			false
		)
	);
	
	dimensionedScalar centerX // x coordinate for center
	(
	   moveFieldsProperties.lookup("centerX")
	);
	
	dimensionedScalar deltaX // mesh size
	(
	   moveFieldsProperties.lookup("deltaX")
	);
	
	dimensionedScalar nMoveCellThreshHold // trigger for move
	(
	   moveFieldsProperties.lookup("nMoveCellThreshHold")
	);
	
	// first of all, we need to find the flame front
	volScalarField magGradT = mag(fvc::grad(T));
	label maxMagGradTi =  findMax(magGradT);
	scalar currentX = mesh.C()[maxMagGradTi].x();
	
	// calculate the move length
	scalar	moveLength = currentX - centerX.value();
	label  moveNCells = std::floor(moveLength/deltaX.value());
	bool flag = false;
	if(mag(moveNCells) > nMoveCellThreshHold.value())
	{
		flag = true;
	}
	bool direaction = false;
	if(mesh.C()[mesh.nCells()].x() > mesh.C()[0].x())
	{
		direaction = true;
	}
	
	// move fields p, he, Yi, U
	if(flag)
	{
		volScalarField& he = thermo.he();
		scalarField pTmp = p.internalField();
		scalarField heTmp = he.internalField();
		vectorField UTmp = U.internalField();
		List<scalarField> YTmp(Y.size());
		forAll(YTmp,i)
		{
			YTmp[i] = Y[i].internalField();
		}
		forAll(mesh.C(),celli)
		{
			label icell;
			if(direaction)
			{
				icell = celli-moveNCells;
			}
			else
			{
				icell = celli+moveNCells;				
			}
			
			if(icell > mesh.C().size()-1)
			{
				icell = mesh.C().size()-1;
			}
			else if(icell < 0)
			{
				icell = 0;
			}
			pTmp[celli] = p[icell];
			heTmp[celli] = he[icell];
			UTmp[celli] = U[icell];
			forAll(YTmp,i)
			{
				YTmp[i][celli] = Y[i][icell]; 
			}
		}
		// update
		forAll(mesh.C(),celli)
		{
			p[celli] = pTmp[celli];
			he[celli] = heTmp[celli];
			U[celli] = UTmp[celli];
			forAll(Y,i)
			{
				Y[i][celli] = YTmp[i][celli];
			}
		}
		// update T, rho, phi and dpdt
		thermo.correct();
		rho = thermo.rho();
		phi = linearInterpolate(rho*U) & mesh.Sf();
		if (thermo.dpdt())
		{
			dpdt = fvc::ddt(p);
		}	
		Info << "%%%%%%%%%%%%%%%%%ATTENTION%%%%%%%%%%%%% "  << endl;
		Info << " move flame front " << endl;
		Info << "%%%%%%%%%%%%%%%%%ATTENTION%%%%%%%%%%%%% "  << endl;
	}		
}