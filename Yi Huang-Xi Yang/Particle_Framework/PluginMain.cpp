#include <maya/MPxCommand.h>
#include <maya/MFnPlugin.h>
#include <maya/MIOStream.h>
#include <maya/MString.h>
#include <maya/MArgList.h>
#include <maya/MGlobal.h>
#include <maya/MSimple.h>
#include <maya/MDoubleArray.h>
#include <maya/MPoint.h>
#include <maya/MPointArray.h>
#include <maya/MFnNurbsCurve.h>
#include <maya/MDGModifier.h>
#include <maya/MPlugArray.h>
#include <maya/MVector.h>
#include <maya/MFnDependencyNode.h>
#include <maya/MStringArray.h>
#include <list>

//#include "LSystemCmd.h"
#include "renderCmd.h"
#include "ParticleNode.h"
//#include "RenderNode.h"
#include <maya/MFnPlugin.h>

MStatus initializePlugin( MObject obj )
{
    MStatus   status = MStatus::kSuccess;
	//MStatus statusN = MStatus::kSuccess;
    MFnPlugin plugin( obj, "MyPlugin", "1.0", "Any");

    // Register Command
    status = plugin.registerCommand( "renderCmd", RenderCmd::creator ,RenderCmd::newSyntax);
	status=plugin.registerNode("ParticleNode", ParticleNode::id, ParticleNode::creator, ParticleNode::initialize);
	//status = plugin.registerNode("RenderNode", RenderNode::id, RenderNode::creator, RenderNode::initialize);
	char buffer[2048];
	sprintf_s(buffer, 2048, "source \"%s/menuInitialization5.mel\";", plugin.loadPath().asChar());
	MGlobal::executeCommand(buffer, true);
    if (!status) {
        status.perror("registerNode");
        return status;
    }
	

    return status;
}

MStatus uninitializePlugin( MObject obj)
{
    MStatus   status = MStatus::kSuccess;
    MFnPlugin plugin( obj );

	status = plugin.deregisterCommand("renderCmd");
	if (!status) {
		status.perror("deregisterCommand");
		return status;
	}

	status= plugin.deregisterNode(ParticleNode::id);
	//status = plugin.deregisterNode(RenderNode::id);
    //delete menu
	MGlobal::executeCommand("if (`window - ex $window`) deleteUI $window;");
	MGlobal::executeCommand("if (`window - ex $renderWindow`) deleteUI $renderWindow;");
	MGlobal::executeCommand("if(`menu -ex MyMenu`) deleteUI MyMenu;");

	if (!status) {
		status.perror("deregisterNode");
		return status;
	}

    return status;
}


