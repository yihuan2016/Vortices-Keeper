#ifndef CreateRenderCmd_H_
#define CreateRenderCmd_H_

#include <maya/MPxCommand.h>
#include <string>
#define MNoVersionString 
#define MNoPluginEntry
#include <maya/MFnPlugin.h>
//#include <maya/MSimple.h>
#include <maya/MGlobal.h>
#include <maya/MArgList.h>
#include <maya/MSyntax.h>
#include <maya/MArgDatabase.h>


class RenderCmd : public MPxCommand
{
public:
	RenderCmd();
	virtual ~RenderCmd();
	static void* creator() { return new RenderCmd(); }
	MStatus doIt(const MArgList& args);
	static MSyntax newSyntax();
	//MStatus parseArgs(const MArgList& args, double &a, double &s, int &ite, MString &gra);
	MStatus parseArgs(const MArgList& args, double &psize, double &density, double &trans,double &shadow,int &stf,int &endf, MString &gra,MString &old);
};



#endif
