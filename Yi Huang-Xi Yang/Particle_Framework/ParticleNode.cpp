#include "ParticleNode.h"
#include <maya/MMatrix.h>
#include <math.h>
#include <string>
#include <vector>
#include <map>
#include "vec.h"
#include <maya/MGlobal.h>
#include <cmath>
#include "tbb/tbb.h"
#include <iostream>

MStatus returnStatus;
MStatus MPxNode::dependsOn(const MPlug &, const MPlug &, bool &) const
{
	return MS::kUnknownParameter;
}

MObject ParticleNode::mTime;
MObject ParticleNode::outputPP;
MTypeId ParticleNode::id(0x80000);
MObject ParticleNode::mCellNumberX;
MObject ParticleNode::mCellNumberY;
MObject ParticleNode::mCellNumberZ;
MObject ParticleNode::mCellSize;
MObject ParticleNode::mTimeStep;
MObject ParticleNode::mTemperature;
MObject ParticleNode::mSmokeDensity;
MObject ParticleNode::mAlpha;
MObject ParticleNode::mBeta;
SmokeSolver3D ParticleNode::g_smokeSolver;
paraSet ParticleNode::paraSets;
std::vector<vector<Vec4f>> ParticleNode::mTracer_stack;
MObject ParticleNode::mFilepath;

//added
MObject ParticleNode::mAdvectionType;
MObject ParticleNode::mEmitterPosX;
MObject ParticleNode::mEmitterPosY;
MObject ParticleNode::mEmitterPosZ;
MObject ParticleNode::mVortConfine_str;
MObject ParticleNode::mEndframe;

MObject ParticleNode::mWindX;
MObject ParticleNode::mWindY;
MObject ParticleNode::mWindZ;
MObject ParticleNode::mParticleNum;

void* ParticleNode::creator()
{
	return new ParticleNode;
}

MStatus ParticleNode::initialize()
{
	MFnUnitAttribute unitAttr;
	MFnTypedAttribute typedAttr;
	MFnNumericAttribute nAttr;
	MStatus returnStatus;

	ParticleNode::mTime = unitAttr.create("time", "tm",
		MFnUnitAttribute::kTime,
		1.0, &returnStatus);
	McheckErr(returnStatus, "ERROR creating ParticleNode time attribute\n");
	//unitAttr.setDefault(1.0);
	MTime maxtime;
	maxtime.setValue(60);
	unitAttr.setMax(maxtime);

	ParticleNode::outputPP = typedAttr.create("outputPosition", "outpp",
		MFnArrayAttrsData::kDynArrayAttrs,
		&returnStatus);
	McheckErr(returnStatus, "ERROR creating ParticleNode output attribute\n");
	//typedAttr.setStorable(false);
	//typedAttr.setKeyable(false);
	//typedAttr.setReadable(true);
	//typedAttr.setWritable(false);
	//**********************output path*************************
	ParticleNode::mFilepath = typedAttr.create("filepath", "fp", MFnData::kString, MFnStringData().create(paraSets.pFile_path), &returnStatus);
	McheckErr(returnStatus, "ERROR creating file path attribute\n");

	//**********************advection control*******************
	ParticleNode::mEndframe = nAttr.create("endframe", "edf", MFnNumericData::kInt, paraSets.pEndframe, &returnStatus);
	nAttr.setMin(1);
	McheckErr(returnStatus, "ERROR creating ParticleNode end frame attribute\n");
	ParticleNode::mAdvectionType = nAttr.create("AdvectionType", "adtype", MFnNumericData::kInt, paraSets.pAdvection_type, &returnStatus);
	nAttr.setMin(-1);
	nAttr.setMax(8);
	McheckErr(returnStatus, "ERROR creating ParticleNode advection type attribute\n");
	ParticleNode::mEmitterPosX = nAttr.create("emitterposX", "emposx", MFnNumericData::kFloat, paraSets.pEmitterPosX, &returnStatus);
	McheckErr(returnStatus, "ERROR creating emitter x attribute\n");
	ParticleNode::mEmitterPosY = nAttr.create("emitterposY", "emposy", MFnNumericData::kFloat, paraSets.pEmitterPosY, &returnStatus);
	McheckErr(returnStatus, "ERROR creating emitter y attribute\n");
	ParticleNode::mEmitterPosZ = nAttr.create("emitterposZ", "emposz", MFnNumericData::kFloat, paraSets.pEmitterPosZ, &returnStatus);
	McheckErr(returnStatus, "ERROR creating emitter z attribute\n");
	ParticleNode::mVortConfine_str = nAttr.create("vortconfine", "vortconf", MFnNumericData::kFloat, paraSets.pVort_confine_str, &returnStatus);
	McheckErr(returnStatus, "ERROR creating vort confine attribute\n");
	//*********************wind control*************************
	ParticleNode::mWindX = nAttr.create("WindX", "wx", MFnNumericData::kFloat, paraSets.pWindX, &returnStatus);
	McheckErr(returnStatus, "ERROR creating wind x attribute\n");
	ParticleNode::mWindY = nAttr.create("WindY", "wy", MFnNumericData::kFloat, paraSets.pWindY, &returnStatus);
	McheckErr(returnStatus, "ERROR creating wind y attribute\n");
	ParticleNode::mWindZ = nAttr.create("WindZ", "wz", MFnNumericData::kFloat, paraSets.pWindZ, &returnStatus);
	McheckErr(returnStatus, "ERROR creating wind z attribute\n");
	ParticleNode::mParticleNum = nAttr.create("particleNumber", "pnum", MFnNumericData::kInt, paraSets.pParticleNum, &returnStatus);
	McheckErr(returnStatus, "ERROR creating particle number attribute\n");

	//**********************cell number*******************
	ParticleNode::mCellNumberX = nAttr.create("cellNumberX", "cNumX", MFnNumericData::kInt, paraSets.pCellNumberX, &returnStatus);
	nAttr.setMin(0);
	//nAttr.setMax(360.0);
	McheckErr(returnStatus, "ERROR creating ParticleNode cell number x attribute\n");
	ParticleNode::mCellNumberY = nAttr.create("cellNumberY", "cNumY", MFnNumericData::kInt, paraSets.pCellNumberY, &returnStatus);
	nAttr.setMin(0);
	//nAttr.setMax(360.0);
	McheckErr(returnStatus, "ERROR creating ParticleNode cell number y attribute\n");
	ParticleNode::mCellNumberZ = nAttr.create("cellNumberZ", "cNumZ", MFnNumericData::kInt, paraSets.pCellNumberZ, &returnStatus);
	nAttr.setMin(0);
	//nAttr.setMax(360.0);
	McheckErr(returnStatus, "ERROR creating ParticleNode cell number z attribute\n");

	//*********************size of cell********************
	ParticleNode::mCellSize = nAttr.create("cellSize", "cSize", MFnNumericData::kFloat, paraSets.pCellSize, &returnStatus);
	nAttr.setMin(0.0);
	//nAttr.setMax(50.0);
	//nAttr.setDefault(20.0);
	McheckErr(returnStatus, "ERROR creating ParticleNode cell size attribute\n");

	//*********************time step************************
	ParticleNode::mTimeStep = nAttr.create("timeStep", "tStep", MFnNumericData::kFloat, paraSets.pTimeStep, &returnStatus);
	McheckErr(returnStatus, "ERROR creating ParticleNode time step attribute\n");

	//********************emitter temperature****************
	ParticleNode::mTemperature = nAttr.create("emitterTemperature", "eTemp", MFnNumericData::kFloat, paraSets.pTemperature, &returnStatus);
	McheckErr(returnStatus, "ERROR creating ParticleNode grammar attribute\n");

	//********************smoke control**********************
	ParticleNode::mSmokeDensity = nAttr.create("smokeDensity", "dens", MFnNumericData::kFloat, paraSets.pSmokeDensity, &returnStatus);
	McheckErr(returnStatus, "ERROR creating ParticleNode smoke density attribute\n");

	ParticleNode::mAlpha = nAttr.create("Alpha", "alpha", MFnNumericData::kFloat, paraSets.pAlpha, &returnStatus);
	McheckErr(returnStatus, "ERROR creating ParticleNode alpha attribute\n");

	ParticleNode::mBeta = nAttr.create("Beta", "beta", MFnNumericData::kFloat, paraSets.pBeta, &returnStatus);
	McheckErr(returnStatus, "ERROR creating ParticleNode beta attribute\n");

	//*************************add attributes*********************

	returnStatus = addAttribute(ParticleNode::mTime);
	McheckErr(returnStatus, "ERROR adding time attribute\n");
	returnStatus = addAttribute(ParticleNode::mEndframe);
	returnStatus = addAttribute(ParticleNode::outputPP);
	McheckErr(returnStatus, "ERROR adding outputPP attribute\n");
	returnStatus = addAttribute(ParticleNode::mCellNumberX);
	McheckErr(returnStatus, "ERROR adding Cell Number X attribute\n");
	returnStatus = addAttribute(ParticleNode::mCellNumberY);
	McheckErr(returnStatus, "ERROR adding Cell Number Y attribute\n");
	returnStatus = addAttribute(ParticleNode::mCellNumberZ);
	McheckErr(returnStatus, "ERROR adding Cell Number Z attribute\n");
	returnStatus = addAttribute(ParticleNode::mCellSize);
	McheckErr(returnStatus, "ERROR adding cell size attribute\n");
	returnStatus = addAttribute(ParticleNode::mTimeStep);
	McheckErr(returnStatus, "ERROR adding time step attribute\n");
	returnStatus = addAttribute(ParticleNode::mTemperature);
	McheckErr(returnStatus, "ERROR adding temperature attribute\n");
	returnStatus = addAttribute(ParticleNode::mSmokeDensity);
	McheckErr(returnStatus, "ERROR adding smoke density attribute\n");
	returnStatus = addAttribute(ParticleNode::mAlpha);
	McheckErr(returnStatus, "ERROR adding alpha attribute\n");
	returnStatus = addAttribute(ParticleNode::mBeta);
	McheckErr(returnStatus, "ERROR adding beta attribute\n");

	//******************added********************************
	returnStatus = addAttribute(ParticleNode::mAdvectionType);
	returnStatus = addAttribute(ParticleNode::mVortConfine_str);
	returnStatus = addAttribute(ParticleNode::mEmitterPosX);
	returnStatus = addAttribute(ParticleNode::mEmitterPosY);
	returnStatus = addAttribute(ParticleNode::mEmitterPosZ);
	McheckErr(returnStatus, "ERROR adding advection attribute\n");
	returnStatus = addAttribute(ParticleNode::mFilepath);
	McheckErr(returnStatus, "ERROR adding file path attribute\n");
	returnStatus = addAttribute(ParticleNode::mWindX);
	returnStatus = addAttribute(ParticleNode::mWindY);
	returnStatus = addAttribute(ParticleNode::mWindZ);
	McheckErr(returnStatus, "ERROR adding wind attribute\n");
	returnStatus = addAttribute(ParticleNode::mParticleNum);
	McheckErr(returnStatus, "ERROR adding particle number attribute\n");
	//**********************add affects***********************
	returnStatus = attributeAffects(ParticleNode::mTime,
		ParticleNode::outputPP);
	McheckErr(returnStatus, "ERROR in attributeAffects\n");
	returnStatus = attributeAffects(ParticleNode::mTimeStep,
		ParticleNode::outputPP);
	McheckErr(returnStatus, "ERROR in attributeAffects\n");
	returnStatus = attributeAffects(ParticleNode::mCellNumberX,
		ParticleNode::outputPP);
	McheckErr(returnStatus, "ERROR in attributeAffects\n");
	returnStatus = attributeAffects(ParticleNode::mCellNumberY,
		ParticleNode::outputPP);
	McheckErr(returnStatus, "ERROR in attributeAffects\n");
	returnStatus = attributeAffects(ParticleNode::mCellNumberZ,
		ParticleNode::outputPP);
	McheckErr(returnStatus, "ERROR in attributeAffects\n");
	returnStatus = attributeAffects(ParticleNode::mCellSize,
		ParticleNode::outputPP);
	McheckErr(returnStatus, "ERROR in attributeAffects\n");
	returnStatus = attributeAffects(ParticleNode::mTemperature,
		ParticleNode::outputPP);
	McheckErr(returnStatus, "ERROR in attributeAffects\n");
	returnStatus = attributeAffects(ParticleNode::mSmokeDensity,
		ParticleNode::outputPP);
	McheckErr(returnStatus, "ERROR in attributeAffects\n");
	returnStatus = attributeAffects(ParticleNode::mAlpha,
		ParticleNode::outputPP);
	McheckErr(returnStatus, "ERROR in attributeAffects\n");
	returnStatus = attributeAffects(ParticleNode::mBeta,
		ParticleNode::outputPP);
	McheckErr(returnStatus, "ERROR in attributeAffects\n");

	//**************added*****************************
	returnStatus = attributeAffects(ParticleNode::mAdvectionType,
		ParticleNode::outputPP);
	McheckErr(returnStatus, "ERROR in attributeAffects\n");

	returnStatus = attributeAffects(ParticleNode::mEmitterPosX,
		ParticleNode::outputPP);
	McheckErr(returnStatus, "ERROR in attributeAffects\n");
	returnStatus = attributeAffects(ParticleNode::mEmitterPosY,
		ParticleNode::outputPP);
	McheckErr(returnStatus, "ERROR in attributeAffects\n");

	returnStatus = attributeAffects(ParticleNode::mEmitterPosZ,
		ParticleNode::outputPP);
	McheckErr(returnStatus, "ERROR in attributeAffects\n");

	returnStatus = attributeAffects(ParticleNode::mVortConfine_str,
		ParticleNode::outputPP);
	McheckErr(returnStatus, "ERROR in attributeAffects\n");

	returnStatus = attributeAffects(ParticleNode::mFilepath,
		ParticleNode::outputPP);
	McheckErr(returnStatus, "ERROR in attributeAffects\n");
	returnStatus = attributeAffects(ParticleNode::mEndframe,
		ParticleNode::outputPP);
	McheckErr(returnStatus, "ERROR in attributeAffects\n");

	returnStatus = attributeAffects(ParticleNode::mWindX,
		ParticleNode::outputPP);
	McheckErr(returnStatus, "ERROR in attributeAffects\n");
	returnStatus = attributeAffects(ParticleNode::mWindY,
		ParticleNode::outputPP);
	McheckErr(returnStatus, "ERROR in attributeAffects\n");
	returnStatus = attributeAffects(ParticleNode::mWindZ,
		ParticleNode::outputPP);
	McheckErr(returnStatus, "ERROR in attributeAffects\n");
	returnStatus = attributeAffects(ParticleNode::mParticleNum,
		ParticleNode::outputPP);
	McheckErr(returnStatus, "ERROR in attributeAffects\n");

	//initialize smoke solver
	//returnStatus = updateSolver();
	//McheckErr(returnStatus, "ERROR initializing smoke solver\n");
	//returnStatus = ComputeTracerStack();
	//McheckErr(returnStatus, "ERROR initializing smoke solver\n");

	return MS::kSuccess;
}


MStatus ParticleNode::compute(const MPlug& plug, MDataBlock& data)

{
	MStatus returnStatus;


	if (plug == outputPP) {
		/* Get time */
		MDataHandle timeData = data.inputValue(mTime, &returnStatus);
		McheckErr(returnStatus, "Error getting time data handle\n");
		MTime nTime = timeData.asTime();
		//****************get input value*****************************
		MDataHandle timeStepData = data.inputValue(mTimeStep, &returnStatus);
		McheckErr(returnStatus, "Error getting time step data handle\n");
		float nTimeStep = timeStepData.asFloat();

		MDataHandle numXData = data.inputValue(mCellNumberX, &returnStatus);
		McheckErr(returnStatus, "Error getting cell number x data handle\n");
		int nCellNumberX = numXData.asInt();

		MDataHandle numYData = data.inputValue(mCellNumberY, &returnStatus);
		McheckErr(returnStatus, "Error getting cell number x data handle\n");
		int nCellNumberY = numYData.asInt();

		MDataHandle numZData = data.inputValue(mCellNumberZ, &returnStatus);
		McheckErr(returnStatus, "Error getting cell number x data handle\n");
		int nCellNumberZ = numZData.asInt();

		MDataHandle tempData = data.inputValue(mTemperature, &returnStatus);
		McheckErr(returnStatus, "Error getting temperature data handle\n");
		float nTemperature = tempData.asFloat();

		MDataHandle densityData = data.inputValue(mSmokeDensity, &returnStatus);
		McheckErr(returnStatus, "Error getting smoke density data handle\n");
		float nSmokeDensity = densityData.asFloat();

		MDataHandle sizeData = data.inputValue(mCellSize, &returnStatus);
		McheckErr(returnStatus, "Error getting smoke cell size data handle\n");
		float nCellSize = sizeData.asFloat();

		MDataHandle alphaData = data.inputValue(mAlpha, &returnStatus);
		McheckErr(returnStatus, "Error getting smoke alpha data handle\n");
		float nAlpha = alphaData.asFloat();

		MDataHandle betaData = data.inputValue(mBeta, &returnStatus);
		McheckErr(returnStatus, "Error getting smoke beta data handle\n");
		float nBeta = betaData.asFloat();

		//advection control
		MDataHandle adtypedata = data.inputValue(mAdvectionType, &returnStatus);
		McheckErr(returnStatus, "Error getting advection type data handle\n");
		int nAdvectionType = adtypedata.asInt();
		MDataHandle emitterxdata = data.inputValue(mEmitterPosX, &returnStatus);
		McheckErr(returnStatus, "Error getting emitter x data handle\n");
		float nEmitterPosX = emitterxdata.asFloat();
		MDataHandle emitterydata = data.inputValue(mEmitterPosY, &returnStatus);
		McheckErr(returnStatus, "Error getting emitter y data handle\n");
		float nEmitterPosY = emitterydata.asFloat();
		MDataHandle emitterzdata = data.inputValue(mEmitterPosZ, &returnStatus);
		McheckErr(returnStatus, "Error getting emitter z data handle\n");
		float nEmitterPosZ = emitterzdata.asFloat();
		MDataHandle vortdata = data.inputValue(mVortConfine_str, &returnStatus);
		McheckErr(returnStatus, "Error getting emitter x data handle\n");
		float nVortConfine_str = vortdata.asFloat();

		MDataHandle edData = data.inputValue(mEndframe, &returnStatus);
		McheckErr(returnStatus, "Error getting end frame data handle\n");
		int nEndframe = edData.asInt();

		//output file path
		MDataHandle pathData = data.inputValue(mFilepath, &returnStatus);
		McheckErr(returnStatus, "Error getting file path data handle\n");
		MString nFilepath = pathData.asString();

		//wind
		MDataHandle particlenumdata = data.inputValue(mParticleNum, &returnStatus);
		McheckErr(returnStatus, "Error getting particle num data handle\n");
		int nParticleNum = particlenumdata.asInt();
		MDataHandle windxdata = data.inputValue(mWindX, &returnStatus);
		McheckErr(returnStatus, "Error getting wind x data handle\n");
		float nWindX = windxdata.asFloat();
		MDataHandle windydata = data.inputValue(mWindY, &returnStatus);
		McheckErr(returnStatus, "Error getting wind y data handle\n");
		float nWindY = windydata.asFloat();
		MDataHandle windzdata = data.inputValue(mWindZ, &returnStatus);
		McheckErr(returnStatus, "Error getting wind z data handle\n");
		float nWindZ = windzdata.asFloat();

		//****************get output attribute*******************************

		MDataHandle outputHandle = data.outputValue(outputPP, &returnStatus);
		McheckErr(returnStatus, "ERROR getting output data handle\n");
		MFnArrayAttrsData dataCreator;
		MObject newOutputData = dataCreator.create(&returnStatus);
		McheckErr(returnStatus, "ERROR creating outputData");
		MVectorArray positionArray = dataCreator.vectorArray("position", &returnStatus);
		McheckErr(returnStatus, "ERROR getting position array\n");
		MDoubleArray idArray = dataCreator.doubleArray("id", &returnStatus);
		McheckErr(returnStatus, "ERROR getting id array\n");

		MVectorArray scaleArray = dataCreator.vectorArray("scale", &returnStatus);
		McheckErr(returnStatus, "ERROR getting scale array\n");


		//check solver parameters change
		//if changed, start over
		bool change = false;
		bool recomputesolver = false;
		if (nCellNumberX - paraSets.pCellNumberX != 0)
		{
			paraSets.pCellNumberX = nCellNumberX;
			change = true;
			if (recomputesolver == false) recomputesolver = true;
			
		} 
		if (nCellNumberY - paraSets.pCellNumberY != 0)
		{
			paraSets.pCellNumberY = nCellNumberY;
			change = true;
			if (recomputesolver == false) recomputesolver = true;
			
		}
		if (nCellNumberZ - paraSets.pCellNumberZ != 0)
		{
			paraSets.pCellNumberZ = nCellNumberZ;
			change = true;
			if (recomputesolver == false) recomputesolver = true;
		}
		if (fabs(nTimeStep-paraSets.pTimeStep)>FLT_EPSILON)
		{
			paraSets.pTimeStep = nTimeStep;
			change = true;
			if (recomputesolver == false) recomputesolver = true;
		}
		if (fabs(nTemperature-paraSets.pTemperature)>FLT_EPSILON)
		{
			paraSets.pTemperature = nTemperature;
			change = true;
			if (recomputesolver == false) recomputesolver = true;
		}
		if (fabs(nSmokeDensity - paraSets.pSmokeDensity)>FLT_EPSILON)
		{
			paraSets.pSmokeDensity = nSmokeDensity;
			change = true;
			if (recomputesolver == false) recomputesolver = true;
		}
		if (fabs(nCellSize - paraSets.pCellSize)>FLT_EPSILON)
		{
			paraSets.pCellSize = nCellSize;
			change = true;
			if (recomputesolver == false) recomputesolver = true;
		}
		if (fabs(nAlpha - paraSets.pAlpha)>FLT_EPSILON)
		{
			paraSets.pAlpha = nAlpha;
			change = true;
			if (recomputesolver == false) recomputesolver = true;
		}
		if (fabs(nBeta - paraSets.pBeta)>FLT_EPSILON)
		{
			paraSets.pBeta = nBeta;
			change = true;
			if (recomputesolver == false) recomputesolver = true;
		}
		if (fabs(nEmitterPosX - paraSets.pEmitterPosX)>FLT_EPSILON)
		{
			paraSets.pEmitterPosX = nEmitterPosX;
			change = true;
		}
		if (fabs(nEmitterPosY - paraSets.pEmitterPosY)>FLT_EPSILON)
		{
			paraSets.pEmitterPosY = nEmitterPosY;
			change = true;
		}
		if (fabs(nEmitterPosZ - paraSets.pEmitterPosZ)>FLT_EPSILON)
		{
			paraSets.pEmitterPosZ = nEmitterPosZ;
			change = true;
		}
		if (fabs(nVortConfine_str - paraSets.pVort_confine_str)>FLT_EPSILON)
		{
			paraSets.pVort_confine_str = nVortConfine_str;
			change = true;
			if (recomputesolver == false) recomputesolver = true;
		}
		if (abs(nAdvectionType - paraSets.pAdvection_type)>0)
		{
			paraSets.pAdvection_type = nAdvectionType;
			change = true;
			if (recomputesolver == false) recomputesolver = true;
		}
		if (abs(nEndframe - paraSets.pEndframe)>0)
		{
			paraSets.pEndframe = nEndframe;
			change = true;
			if (recomputesolver == false) recomputesolver = true;
		}
		if (abs(nParticleNum - paraSets.pParticleNum)>0)
		{
			paraSets.pParticleNum = nParticleNum;
			change = true;
			if (recomputesolver == false) recomputesolver = true;
		}
		if (fabs(nWindX - paraSets.pWindX)>FLT_EPSILON)
		{
			paraSets.pWindX = nWindX;
			change = true;
			if (recomputesolver == false) recomputesolver = true;
		}
		if (fabs(nWindY - paraSets.pWindY)>FLT_EPSILON)
		{
			paraSets.pWindY = nWindY;
			change = true;
			if (recomputesolver == false) recomputesolver = true;
		}
		if (fabs(nWindZ - paraSets.pWindZ)>FLT_EPSILON)
		{
			paraSets.pWindZ = nWindZ;
			change = true;
			if (recomputesolver == false) recomputesolver = true;
		}
		if(strcmp(paraSets.pFile_path, nFilepath.asChar()) != 0)
		//if (paraSets.pFile_path != nFilepath.asChar())
		{
			
			int n = sprintf(paraSets.pFile_path, "%s", nFilepath.asChar());
			if (recomputesolver == false) returnStatus=ComputeTracerStack();

			
		}
		MAnimControl aniControl;
		if (change)
		{
			
			MTime restartTime;
			restartTime.setValue(0);
			aniControl.setCurrentTime(restartTime);			
		}
		if (recomputesolver == true)
		{
			returnStatus = updateSolver();
			returnStatus = ComputeTracerStack();
		}
		//**************************IVOCK main*****************
		
		/*
		//int frame = (int)nTime.as(MTime::kFilm);
		MTime animTime = aniControl.currentTime();
		//int frame = (int)animTime.as(MTime::kFilm);
		frame++;
		//MGlobal::displayInfo((std::to_string(frame)).c_str());
		float T = 0.1*(float)frame;

		//g_smokeSolver.setSmoke(10.0, 0.15, 0.25, 500, 10.0);
		g_smokeSolver.setSmoke(10.0, paraSets.pAlpha, paraSets.pBeta, paraSets.pTemperature, paraSets.pSmokeDensity);
		g_smokeSolver.set_vort_confine(paraSets.pVort_confine_str);
		g_smokeSolver.time_step(paraSets.pTimeStep, paraSets.pAdvection_type);

		paraSets.pCurrentFrame = frame;
		int tracerNum = g_smokeSolver.tracers.size();
		int n = sprintf(file_path, "%s", "D:/MayaOut");
		std::cout << file_path << endl;
		//MGlobal::displayInfo((std::to_string(tracerNum)).c_str());
		for (int i = 0; i < tracerNum; i++)
		{

			MVector pos;
			pos[0] = g_smokeSolver.tracers[i][0];
			pos[1] = g_smokeSolver.tracers[i][1];
			pos[2] = g_smokeSolver.tracers[i][2];
			MVector scale;
			scale[0] = 0.01;
			scale[1] = 0.01;
			scale[2] = 0.01;
			positionArray.append(pos);
			idArray.append((double)i);
			scaleArray.append(scale);

			
		}
		//////////////////////////////
		//add by Yi Huang
		g_smokeSolver.time_step(paraSets.pTimeStep, paraSets.pAdvection_type);
		g_smokeSolver.write_tracers(file_path, frame);
		printf("frame %d done\n", frame);

		*/
		//***************change update way****by Xi Yang**********
		MTime animTime = aniControl.currentTime();
		int frame = (int)animTime.as(MTime::kFilm);
		if (frame >= paraSets.pEndframe)
		{
			MTime restartTime;
			restartTime.setValue(0);
			aniControl.setCurrentTime(restartTime);
			frame = 0;
		}
		int tracerNum = mTracer_stack[frame].size();
		for (int i = 0; i < tracerNum; i++)
		{

			MVector pos;
			pos[0] = mTracer_stack[frame][i][0]+paraSets.pEmitterPosX;
			pos[1] = mTracer_stack[frame][i][1]+paraSets.pEmitterPosY;
			pos[2] = mTracer_stack[frame][i][2]+paraSets.pEmitterPosZ;
			MVector scale;
			scale[0] = 0.1;
			scale[1] = 0.1;
			scale[2] = 0.1;
			positionArray.append(pos);
			idArray.append((double)i);
			scaleArray.append(scale);


		}
		//*****************************************************

		//MGlobal::displayInfo("computation completed");
		outputHandle.set(newOutputData);
		//MGlobal::displayInfo("set output");
		data.setClean(plug);


		
	}
	else
		return MS::kUnknownParameter;

	return MS::kSuccess;
}

MStatus ParticleNode::updateSolver()
{
	
	//************************initialize ivock solver**********
	g_smokeSolver.tracers.resize(0);
	paraSets.pCurrentFrame = 0;
	//advection_type = paraSets.pAdvection_type;  //0~8
	//vort_confine_str = paraSets.pVort_confine_str;

	


	int nx = paraSets.pCellNumberX, ny = paraSets.pCellNumberY, nz = paraSets.pCellNumberZ;
	float L = paraSets.pCellSize;
	float g_h = L / (float)nx;
	g_smokeSolver.init(nx, ny, nz, L);
	g_smokeSolver.setSmoke(10.0, paraSets.pAlpha, paraSets.pBeta, paraSets.pTemperature, paraSets.pSmokeDensity);
	

	buffer3Dc bc_des;
	bc_des.init(nx, ny, nz);
	bc_des.setZero();
	for (int k = 0; k<nz; k++)for (int j = 0; j<ny; j++)for (int i = 0; i<nx; i++)
	{
		//0:fluid;1:air;2:solid
		if (i<1) bc_des(i, j, k) = 1;
		if (j<1) bc_des(i, j, k) = 2;
		if (k<1) bc_des(i, j, k) = 1;

		if (i >= nx - 1) bc_des(i, j, k) = 1;
		if (j >= ny - 1) bc_des(i, j, k) = 1;
		if (k >= nz - 1) bc_des(i, j, k) = 1;

		float x = g_h * i;
		float y = g_h * j;
		float z = g_h * k;
		float X = g_h * (float)nx;
		float Y = g_h * (float)ny;
		float Z = g_h * (float)nz;
		//if (sqrt((x-0.5*X)*(x-0.5*X)+(y-0.5*Y)*(y-0.5*Y)+(z-0.5*Z)*(z-0.5*Z))<=0.1*L)
		if (sqrt((x - 0.52*L)*(x - 0.52*L) + (y - 0.125*L)*(y - 0.125*L) + (z - 0.5*L)*(z - 0.5*L)) <= 0.02*L)
		{
			bc_des(i, j, k) = 2;
		}

	}
	g_smokeSolver.set_boundary(bc_des);

	bc_des.setZero();
	for (int k = 0; k<nz; k++)for (int j = 0; j<ny; j++)for (int i = 0; i<nx; i++)
	{
		float x = g_h * i;
		float y = g_h * j;
		float z = g_h * k;
		//1:source; 2:clear
		if (sqrt((x - 0.5*L)*(x - 0.5*L) + (y - 0.12*L)*(y - 0.12*L) + (z - 0.5*L)*(z - 0.5*L)) <= 0.07*L)
		{
			bc_des(i, j, k) = 1;
		}
	}
	g_smokeSolver.set_heat(bc_des);
	//g_smokeSolver.setEmitter(Vec3f(0.5*L, 0.12*L, 0.5*L), 0.015*L, 12288);
	g_smokeSolver.setEmitter(Vec3f(0.5*L, 0.12*L, 0.5*L), 0.015*L, paraSets.pParticleNum);
	//********************************************************
	
	return MS::kSuccess;
}

MStatus ParticleNode::ComputeTracerStack()
{
	//MGlobal::displayInfo((std::to_string(paraSets.pParticleNum)).c_str());
	for (int i = 0; i < mTracer_stack.size(); i++)
	{
		mTracer_stack[i].clear();
	}
	mTracer_stack.clear();
	paraSets.pWindPos.clear();
	//initial pos
	vec3 WindAcc(paraSets.pWindX, paraSets.pWindY, paraSets.pWindZ);
	bool useWind = true;
	if (fabs(WindAcc[0]) < FLT_EPSILON && fabs(WindAcc[1]) < FLT_EPSILON && fabs(WindAcc[2]) < FLT_EPSILON)
		useWind = false;
	if (useWind) {
		for (int frame = paraSets.pStartframe; frame < paraSets.pEndframe; frame++)
		{
			vec3 formerpos(0, 0, 0);
			if (frame > paraSets.pStartframe)
				formerpos = paraSets.pWindPos[frame - 1];
			vec3 pos = WindAcc*frame*paraSets.pTimeStep*paraSets.pTimeStep + formerpos;
			paraSets.pWindPos.push_back(pos);
		}
	}

	for (int frame = paraSets.pStartframe; frame < paraSets.pEndframe; frame++)
	{

		//MGlobal::displayInfo((std::to_string(frame)).c_str());
		float T = 0.1*(float)frame;

		//g_smokeSolver.setSmoke(10.0, 0.15, 0.25, 500, 10.0);
		g_smokeSolver.setSmoke(10.0, paraSets.pAlpha, paraSets.pBeta, paraSets.pTemperature, paraSets.pSmokeDensity);
		g_smokeSolver.set_vort_confine(paraSets.pVort_confine_str);
		g_smokeSolver.time_step(paraSets.pTimeStep, paraSets.pAdvection_type);


		int tracerNum = g_smokeSolver.tracers.size();
		vector<Vec4f> currentTracer;
		for (int i = 0; i < tracerNum; i++)
		{
			Vec4f pos;
			if (useWind) {
				pos[0] = g_smokeSolver.tracers[i][0] + paraSets.pWindPos[g_smokeSolver.tracer_life[i]][0];
				pos[1] = g_smokeSolver.tracers[i][1] + paraSets.pWindPos[g_smokeSolver.tracer_life[i]][0];
				pos[2] = g_smokeSolver.tracers[i][2] + paraSets.pWindPos[g_smokeSolver.tracer_life[i]][0];
				pos[3] = g_smokeSolver.tracers[i][3];
			}
			else
			{
				pos[0] = g_smokeSolver.tracers[i][0];
				pos[1] = g_smokeSolver.tracers[i][1];
				pos[2] = g_smokeSolver.tracers[i][2];
				pos[3] = g_smokeSolver.tracers[i][3];
			}
			currentTracer.push_back(pos);
		}
		mTracer_stack.push_back(currentTracer);
		

		//int n = sprintf(file_path, "%s", "D:/MayaOut");
		//std::cout << file_path << endl;
		//MGlobal::displayInfo((std::to_string(tracerNum)).c_str());
		//MGlobal::displayInfo((std::to_string(g_smokeSolver.tracer_life.size())).c_str());

		//mTracer_stack.push_back(g_smokeSolver.tracers);
		//////////////////////////////
		//add by Yi Huang
		//g_smokeSolver.write_tracers(file_path, frame);
		if (strcmp(paraSets.pFile_path, "Please select a directory") != 0)
			//g_smokeSolver.write_tracers(paraSets.pFile_path, frame);
			g_smokeSolver.write_tracers(paraSets.pFile_path, frame, currentTracer);
		printf("frame %d done\n", frame);
	}
	return MS::kSuccess;
}

/*createNode LSystemNode - n "Lnode";
polySphere - n "lsphere";
instancer - n "lflower";
connectAttr lsphere.matrix lflower.inputHierarchy[0];
connectAttr time1.outTime Lnode.time;
connectAttr Lnode.outpp lflower.inputPoints;*/