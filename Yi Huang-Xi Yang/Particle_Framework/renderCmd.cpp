
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
#include "renderCmd.h"
#include <cstdio>

#include "particle_renderer\particle_renderer\read_particles.h"
#include "particle_renderer\particle_renderer\render.h"
#include "particle_renderer\particle_renderer\shadow.h"
#include "particle_renderer\particle_renderer\write_sgi.h"
#include "particle_renderer\particle_renderer\obj_reader.h"

#include "particle_renderer\particle_renderer\read_particles.cpp"
#include "particle_renderer\particle_renderer\render.cpp"
#include "particle_renderer\particle_renderer\shadow.cpp"
#include "particle_renderer\particle_renderer\write_sgi.cpp"
#include "particle_renderer\particle_renderer\obj_reader.cpp"
#include "particle_renderer\particle_renderer\bfstream.cpp"

#include <maya/MGlobal.h>
#include <list>

//0.008(particle size) 1(density scale) 0.006(particle transparency)  3 10 -5(light position) 10(shadow map quality) 0(start frame) 100(end frame)



const char *particleSize_flag = "-ps", *particleSize_longflag = "-particleSize";
const char *density_flag = "d", *density_longflag = "-density";
const char *transparancy_flag = "ts", *transparancy_longflag = "transparency";
const char *shadow_flag = "sd", *shadow_longflag = "shadowQuality";
const char *startframe_flag = "sf", *startframe_longflag = "startframe";
const char *endframe_flag = "ef", *endframe_longflag = "endframe";
const char *rendefile_flag = "rdf", *rendfile_longflag = "rendefile";
const char *binfile_flag = "bin", *binfile_longflag = "binfile";

double pi = 3.14159265359;
double radius = 0.08, opacity = 0.06;
Vec3f light_position(10, -3, -3);
double density_scale = 0.7;
int startFrame = 1;
int endFrame = 100;
double quality_factor = 10;

MString file = "D:/MayaOut";
MString binFile = "D:/MayaOut";
double frand2(double a, double b)
{
	return (((double)rand()) / ((double)RAND_MAX)*(b - a)) + a;
}
RenderCmd::RenderCmd() : MPxCommand()
{
}

RenderCmd::~RenderCmd()
{
}

MStatus RenderCmd::doIt(const MArgList& args)
{
	// message in Maya output window
	//cout<<"Implement Me!"<<endl;
	//std::cout.flush();
	//MString test1 = "curve -d 1 -p 0 0 0 -p 10 10 0 -k 0 -k 1 -name curve1;";
	//MString test2 = "select -r nurbsCircle1 curve1;";
	//MString test3 = "extrude -ch true -rn false -po 1 -et 2 -ucp 1 -fpt 1 -upn 1 -rotation 0 -scale 1 -rsp 1 \"nurbsCircle1\" \"curve1\";";
	//MGlobal::executeCommand(test1);
	//MGlobal::executeCommand(test2);
	//MGlobal::executeCommand(test3);

	MStatus status = parseArgs(args, radius, density_scale,opacity, quality_factor,startFrame,endFrame,file,binFile);

	//doing things
	radius=radius/100.0f;
	opacity=opacity/100.0f;
	Array1<Vec3f> x, base_rgb, lit_rgb;
	Array1<float> d, temp;
	Vec3f light_position(10, -3, -3);
	
	
	const char *fileFormat = file.asChar();
	const char *binf = binFile.asChar();
	std::printf("radius %f, opacity %f, light (%f %f %f)\n", radius, opacity,
		light_position[0], light_position[1], light_position[2]);

	Array1<float> illum;
	Array2f shadow_map;
	Array2<Vec3f> image(1440, 960);
	Array1<Vec3f> ball;


	ball.resize(163840);
	int num_get = 0;
	while (num_get<163840)
	{
		float x = frand2(-0.2, 0.2);
		float y = frand2(-0.2, 0.2);
		float z = frand2(-0.2, 0.2);
		if (sqrt(x*x + y*y + z*z) <= 0.105&&sqrt(x*x + y*y + z*z) >= 0.095)
		{
			ball[num_get] = Vec3f(x, y + 0.2, z);
			num_get++;
		}
	}

	for (int f = startFrame; f<endFrame; f++) {

		//Vec3f ball_position = Vec3f(0 + 0.005 * 4.0 - 0.01*4.0*(double)f, 0,0);
		Array1<Vec3f> y;

		if (read_particles(x, d, "%s/Particle_data%04d.bin", binf, f)) {

			char filename[256];
			int n = sprintf(filename, "%s/temp_data%04d.bin", binf, f);

			FILE *temp_data = fopen(filename, "rb");
			if (temp_data != NULL)
			{
				temp.resize(x.size());
				size_t result = fread(&(temp[0]), 1, sizeof(float)*temp.size(), temp_data);
				fclose(temp_data);
			}

			y.resize(x.size());
			base_rgb.resize(x.size());
			lit_rgb.resize(x.size());
			for (int p = 0; p < x.size(); p++)
			{
				//y[p]=x[p];
				y[p][0] = x[p][0];
				y[p][1] = x[p][1] + 0.5;
				y[p][2] = x[p][2];

				base_rgb[p] = Vec3f(0.5, 0.7, 0.9);
				lit_rgb[p] = Vec3f(.9, .9, .99);


				if (temp.size() > 0)
				{
					float heat_index = (temp[p]) / 1000.0f;
					heat_index = max(min(heat_index, 1.0f), 0.0f);

					float r = heat_index;
					float g = heat_index*0.4;
					float b = 0.01;

					base_rgb[p] = Vec3f(r, g, b) * 3 + Vec3f(0.08, 0.1, 0.15);
					lit_rgb[p] = Vec3f(r, g, b) * 2 + Vec3f(0.6, 0.6, 0.6);
				}
			}

			for (unsigned int ii = 0; ii < d.size(); ii++)
			{
				d[ii] *= density_scale;
			}

			if (!compute_shadows(y, radius, opacity, d, light_position, quality_factor, illum, shadow_map)) {
				break;
			}

			write_sgi(shadow_map, false, "%s/shadow%04d.sgi", fileFormat, f);

			//for(int i=0;i<60;i++)
			{

				render_smoke(y, illum, base_rgb, lit_rgb, radius, opacity, d, light_position, Vec3f(0.0, 0.0, 0.0),
					Vec3f(0.0, 1.3, -2.5), Vec3f(0, 1.1, 0), 2, image);
				//render_smoke(y, illum, base_rgb, lit_rgb,radius, opacity, d, light_position,  Vec3f(0.0,0.0,0.0), 
				//	Vec3f(0.0, 8.0,-30.5), Vec3f(0,7,0), 2, image);

				//render_smoke(y, illum, base_rgb, lit_rgb,radius, opacity, d, light_position,  Vec3f(0.0,0.0,0.0), 
				//	Vec3f(0, 3.3,-2.6), Vec3f(0,0.9,0), 2, image);
				write_sgi(image, false, "%s/frame%04d.sgi", fileFormat, f);

			}






		}

	}

	// message in scriptor editor
	MGlobal::displayInfo("LSystemCmd excuted!");
	//if (MS::kSuccess != status)
	//	return status;
	return MStatus::kSuccess;
}
/*
const char *angle_flag = "-a", *angle_longflag = "-angle";
const char *grammar_flag = "-g", *grammar_longflag = "-grammar";
const char *iterationNum_flag = "-i", *iterationNum_longflag = "-iteration";

const char *particleSize_flag = "-ps", *particleSize_longflag = "-particleSize";
const char *density_flag = "d", *density_longflag = "-density";
const char *transparancy_flag = "ts", *transparancy_longflag = "transparency";
const char *shadow_flag = "sd", *shadow_longflag = "shadowQuality";
const char *startframe_flag = "sf", *startframe_longflag = "startframe";
const char *endframe_flag = "ef", *endframe_longflag = "endframe";
const char *rendefile_flag = "rdf", *rendfile_longflag = "rendefile";
*/
MSyntax RenderCmd::newSyntax()
{
	MSyntax syntax;
	
	syntax.addFlag(particleSize_flag, particleSize_longflag, MSyntax::kDouble);
	syntax.addFlag(density_flag, density_longflag, MSyntax::kDouble);
	syntax.addFlag(transparancy_flag, transparancy_longflag, MSyntax::kDouble);
	syntax.addFlag(shadow_flag, shadow_longflag, MSyntax::kLong);
	syntax.addFlag(startframe_flag, startframe_longflag, MSyntax::kLong);
	syntax.addFlag(endframe_flag, rendfile_longflag, MSyntax::kLong);
	syntax.addFlag(rendefile_flag, rendfile_longflag, MSyntax::kString);
	syntax.addFlag(binfile_flag, binfile_longflag, MSyntax::kString);
	return syntax;
}
MStatus RenderCmd::parseArgs(const MArgList& args, double &psize, double &density, double &trans, double &shadow, int &stf, int &endf, MString &gra,MString &old)
{

	MStatus status;
	MArgDatabase argData(syntax(), args);
	if (argData.isFlagSet(particleSize_flag)) {
		double tmp;
		status = argData.getFlagArgument(particleSize_flag, 0, tmp);
		if (!status) {
			status.perror("step size flag parsing failed");
			return status;
		}
		psize = tmp;
		MString test;
		test.set(tmp);
		MGlobal::displayInfo(test);
	}
	if (argData.isFlagSet(density_flag)) {
		double des;
		status = argData.getFlagArgument(density_flag, 0, des);
		if (!status) {
			status.perror("angle flag parsing failed");
			return status;
		}
		density = des;
		MString test;
		test.set(des);
		MGlobal::displayInfo(test);
	}
	if (argData.isFlagSet(transparancy_flag)) {
		double itn;
		status = argData.getFlagArgument(transparancy_flag, 0, itn);
		if (!status) {
			status.perror("step size flag parsing failed");
			return status;
		}
		trans = itn;
		MString test;
		test.set(itn);
		MGlobal::displayInfo(test);
	}
	if (argData.isFlagSet(shadow_flag)) {
		int itn;
		status = argData.getFlagArgument(shadow_flag, 0, itn);
		if (!status) {
			status.perror("step size flag parsing failed");
			return status;
		}
		shadow = itn;
		MString test;
		test.set(itn);
		MGlobal::displayInfo(test);
	}
	if (argData.isFlagSet(startframe_flag)) {
		int itn;
		status = argData.getFlagArgument(startframe_flag, 0, itn);
		if (!status) {
			status.perror("step size flag parsing failed");
			return status;
		}
		stf = itn;
		MString test;
		test.set(itn);
		MGlobal::displayInfo(test);
	}
	if (argData.isFlagSet(endframe_flag)) {
		int itn;
		status = argData.getFlagArgument(endframe_flag, 0, itn);
		if (!status) {
			status.perror("step size flag parsing failed");
			return status;
		}
		endf = itn;
		MString test;
		test.set(itn);
		MGlobal::displayInfo(test);
	}



	if (argData.isFlagSet(rendefile_flag)) {
		MString gr;
		status = argData.getFlagArgument(rendefile_flag, 0, gr);
		if (!status) {
			status.perror("step size flag parsing failed");
			return status;
		}
		gra = gr;
		MGlobal::displayInfo(gr);
	}
	if (argData.isFlagSet(binfile_flag)) {
		MString gr;
		status = argData.getFlagArgument(binfile_flag, 0, gr);
		if (!status) {
			status.perror("step size flag parsing failed");
			return status;
		}
		old = gr;
		MGlobal::displayInfo(gr);
	}
	return MS::kSuccess;
}