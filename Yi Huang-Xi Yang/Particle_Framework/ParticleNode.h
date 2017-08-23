#pragma once
#define MNoVersionString 
#define MNoPluginEntry
#include <maya/MFnPlugin.h>

#include <maya/MTime.h>
#include <maya/MFnMesh.h>
#include <maya/MPoint.h>
#include <maya/MFloatPoint.h>
#include <maya/MFloatPointArray.h>
#include <maya/MIntArray.h>
#include <maya/MDoubleArray.h>
#include <maya/MFnUnitAttribute.h>
#include <maya/MFnTypedAttribute.h>
#include <maya/MFnNumericAttribute.h>
#include <maya/MFnStringData.h>
#include <maya/MFnArrayAttrsData.h>
#include <maya/MFnAttribute.h>

#include <maya/MPxNode.h>
#include <maya/MObject.h>
#include <maya/MPlug.h>
#include <maya/MDataBlock.h>
#include <maya/MFnMeshData.h>

#include <maya/MIOStream.h>

#include <maya/MPointArray.h>
#include <maya/MVector.h>
#include <maya/MVectorArray.h>
#include <maya/MIntArray.h>
#include <maya/MDoubleArray.h>
#include <maya/MArrayDataBuilder.h>
#include <maya/MAnimControl.h>

//ivock
#include "Smoke_solver3D.h"
#include "fluid_buffer3D.h"
#include "array.h"
#include "Multigrid3D.h"

#include "vec.h"


//***************

#define McheckErr(stat,msg)			\
	if ( MS::kSuccess != stat ) {	\
		cerr << msg;				\
		return MS::kFailure;		\
	}
class MArrayDataBuilder;
class paraSet;

class ParticleNode : public MPxNode
{
public:
	ParticleNode() {};
	virtual 		~ParticleNode() {};
	virtual MStatus compute(const MPlug& plug, MDataBlock& data);
	static  void*	creator();
	static  MStatus initialize();

	static MObject	mTime;
	static MObject	outputPP;
	static MTypeId	id;

	//parameters for IVOCK
	//number of grid cells
	static MObject mCellNumberX;
	static MObject mCellNumberY;
	static MObject mCellNumberZ;
	//size of cell
	static MObject mCellSize;
	//time step
	static MObject mTimeStep;
	//emitter temperature
	static MObject mTemperature;
	//smoke control
	static MObject mSmokeDensity;
	static MObject mAlpha;
	static MObject mBeta;
	//smoke solver
	static SmokeSolver3D g_smokeSolver;
	
	//initialize
	static MStatus updateSolver();
	static paraSet paraSets;
	

	
	//update per frame
	static std::vector<vector<Vec4f>> mTracer_stack;
	static MStatus ComputeTracerStack();
	static MObject mFilepath;
	//advection type
	static MObject mAdvectionType;
	static MObject mEmitterPosX;
	static MObject mEmitterPosY;
	static MObject mEmitterPosZ;
	static MObject mVortConfine_str;
	static MObject mEndframe;
	//wind force
	static MObject mWindX;
	static MObject mWindY;
	static MObject mWindZ;
	static MObject mParticleNum;


};

class paraSet 
{
public:
	int pCellNumberX=32;
	int pCellNumberY=64;
	int pCellNumberZ=32;
	float pCellSize=5.0;
	float pTimeStep=0.02;
	float pTemperature= 500;
	float pSmokeDensity= 10.0;
	float pAlpha= 0.15;
	float pBeta= 0.25;
	float pVort_confine_str=0.5;
	int pAdvection_type= -1;
	int pCurrentFrame = 0;
	int pParticleNum = 384;
	//*****emitter*****
	float pEmitterPosX = 0.0;
	float pEmitterPosY = 0.0;
	float pEmitterPosZ = 0.0;
	int tracersize = 0;
	int pStartframe=0;
	int pEndframe=60;
	char pFile_path[256]= "D:/MayaOut2";

	//wind force
	float pWindX = 0.0f;
	float pWindY = 0.0f;
	float pWindZ = 0.0f;
	float pMass = 1.0f;
	std::vector<vec3> pWindPos;
	
	/*advection_type	
	case 0: use Semi-Lagrangian advection
	case 1: use SL_iVOCK advection
	case 2: use BFECC advection
	case 3: use BFECC-IVOCK advection
	case 4: use MacCormack advection
	case 5: use MacCormack-IVOCK advection
	case 6: use FLIP advection
	case 7: use FLIP-IVOCK advection
	case 8: use Best Combination
	*/
};