#include "execute.h"
#include <cstdio>



/*
double frand1(double a, double b)
{
	return (((double)rand()) / ((double)RAND_MAX)*(b - a)) + a;
}
*/

void Execute::parseInfor( const char *fileFormat, int startFrame, int endFrame) {
	Array1<Vec3f> x, base_rgb, lit_rgb;
	Array1<float> d, temp;
	float radius = 0.008, opacity = 0.006;
	Vec3f light_position(10, -3, -3);
	float density_scale = 0.7;

	float quality_factor = 10;
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
		float x = frand1(-0.2, 0.2);
		float y = frand1(-0.2, 0.2);
		float z = frand1(-0.2, 0.2);
		if (sqrt(x*x + y*y + z*z) <= 0.105&&sqrt(x*x + y*y + z*z) >= 0.095)
		{
			ball[num_get] = Vec3f(x, y + 1.0, z);
			num_get++;
		}
	}

	for (int f = startFrame; f<endFrame; f++) {

		//Vec3f ball_position = Vec3f(0 + 0.005 * 4.0 - 0.01*4.0*(double)f, 0,0);
		Array1<Vec3f> y;

		if (read_particles(x, d, "%s/Particle_data%04d.bin", fileFormat, f)) {

			char filename[256];
			int n = sprintf(filename, "%s/temp_data%04d.bin", fileFormat, f);

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
}
