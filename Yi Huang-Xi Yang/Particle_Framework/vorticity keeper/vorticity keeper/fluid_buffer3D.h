#ifndef __buffer_h__
#define __buffer_h__

#include "array.h"
#include <vector>

using namespace std;
template<class T>
class Buffer3D {
public:
	Buffer3D() { _nx = _ny = _nz = _n = 0; _hx = _hy = _hz = 0.0; _data = new array1D<T>; }
	~Buffer3D() {}
	array1D<T> *_data;
	vector<vector<vector<T*>>> index;
	uint _blockN;



	// f should be between 0 and 1
	inline void cubic_interp_weights(T f, T& wneg1, T& w0, T& w1, T& w2)
	{
		T f2(f*f), f3(f2*f);
		wneg1 = -T(1. / 3)*f + T(1. / 2)*f2 - T(1. / 6)*f3;
		w0 = 1 - f2 + T(1. / 2)*(f3 - f);
		w1 = f + T(1. / 2)*(f2 - f3);
		w2 = T(1. / 6)*(f3 - f);
	}

	template<class S>
	inline S cubic_interp(const S& value_neg1, const S& value0,
		const S& value1, const S& value2, T f)
	{
		T wneg1, w0, w1, w2;
		cubic_interp_weights(f, wneg1, w0, w1, w2);
		return wneg1*value_neg1 + w0*value0 + w1*value1 + w2*value2;
	}




	uint _nx, _ny, _nz, _n;
	double _hx, _hy, _hz;
	double _ox, _oy, _oz;
	uint _physical_zstride, _physical_ystride, _physical_nx, _physical_ny, _physical_nz, _physical_n, _blockx, _blocky, _blockz;
	void init(uint nx, uint ny, uint nz, double h, double ox, double oy, double oz)
	{
		_nx = nx; _ny = ny; _nz = nz;
		_n = nx*ny*nz;
		_hx = _hy = _hz = h;
		_ox = ox; _oy = oy; _oz = oz;
		_blockN = 8;

		//data are now going to be stored by blocks
		_blockx = (!(_nx%_blockN)) ? _nx / _blockN : (_nx / _blockN + 1);
		_blocky = (!(_ny%_blockN)) ? _ny / _blockN : (_ny / _blockN + 1);
		_blockz = (!(_nz%_blockN)) ? _nz / _blockN : (_nz / _blockN + 1);
		_physical_nx = _blockx*_blockN;
		_physical_ny = _blocky*_blockN;
		_physical_nz = _blockz*_blockN;
		_physical_n = _physical_nx*_physical_ny*_physical_nz;
		_data->alloc(_physical_n);
		_data->setZero();


		index.clear();
		index.resize(_blockz);
		for (int k = 0; k<_blockz; k++)
		{
			index[k].resize(_blocky);
		}
		for (int k = 0; k<_blockz; k++)for (int j = 0; j<_blocky; j++)
		{
			index[k][j].resize(_blockx);
		}
		for (int k = 0; k<_blockz; k++)for (int j = 0; j<_blocky; j++) for (int i = 0; i<_blockx; i++)
		{
			index[k][j][i] = &(_data->getPtr()[((k*_blocky + j)*_blockx + i)*_blockN*_blockN*_blockN]);
		}



	}

	void init(uint nx, uint ny, uint nz)
	{
		_nx = nx; _ny = ny; _nz = nz;
		_n = nx*ny*nz;
		_hx = _hy = _hz = 0;
		_ox = 0; _oy = 0; _oz = 0;
		_blockN = 8;


		_blockx = (!(_nx%_blockN)) ? _nx / _blockN : (_nx / _blockN + 1);
		_blocky = (!(_ny%_blockN)) ? _ny / _blockN : (_ny / _blockN + 1);
		_blockz = (!(_nz%_blockN)) ? _nz / _blockN : (_nz / _blockN + 1);
		_physical_nx = _blockx*_blockN;
		_physical_ny = _blocky*_blockN;
		_physical_nz = _blockz*_blockN;
		_physical_n = _physical_nx*_physical_ny*_physical_nz;

		_data->alloc(_physical_n);
		_data->setZero();

		index.clear();
		index.resize(_blockz);
		for (int k = 0; k<_blockz; k++)
		{
			index[k].resize(_blocky);
		}
		for (int k = 0; k<_blockz; k++)for (int j = 0; j<_blocky; j++)
		{
			index[k][j].resize(_blockx);
		}
		for (int k = 0; k<_blockz; k++)for (int j = 0; j<_blocky; j++) for (int i = 0; i<_blockx; i++)
		{
			index[k][j][i] = &(_data->getPtr()[((k*_blocky + j)*_blockx + i)*_blockN*_blockN*_blockN]);
		}



	}

	Buffer3D<T> initAll(uint nx, uint ny, uint nz) {
		buffer3Dc bc_des;
			bc_des.init(nx, ny, nz);
			bc_des.setZero();
		float L = 20.0;
		float g_h = L / (float)nx;
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
		return bc_des;
	}

	Buffer3D<T> resetAll(uint nx, uint ny, uint nz) {
		buffer3Dc bc_des;
		bc_des.setZero();
		float L = 20.0;
		float g_h = L / (float)nx;
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

		return bc_des;


	}

	array1D<T> *getArray() { return _data; }
	uint getSize() { return _n; }
	void setZero() { _data->setZero(); }
	void free() { _data->free(); index.clear(); }
	void copy(Buffer3D<T> &b)
	{
		memcpy(_data->getPtr(), (b.getArray())->getPtr(), sizeof(T)*_physical_n);

	}
	const T& operator()(uint i, uint j, uint k) const
	{
		uint I = i >> 3, J = j >> 3, K = k >> 3;
		uint ii = i & 7, jj = j & 7, kk = k & 7;
		uint idx = ((K*_blockx*_blocky + J*_blockx + I) << 9) + (kk << 6) + (jj << 3) + ii;
		return (_data->getPtr()[idx]);
		//return (index[K][J][I])[(((kk<<3)+jj)<<3)+ii];
	}
	T& operator()(uint i, uint j, uint k)
	{

		uint I = i >> 3, J = j >> 3, K = k >> 3;
		uint ii = i & 7, jj = j & 7, kk = k & 7;
		uint idx = ((K*_blockx*_blocky + J*_blockx + I) << 9) + (kk << 6) + (jj << 3) + ii;
		return (_data->getPtr()[idx]);
		//return (index[K][J][I])[(((kk<<3)+jj)<<3)+ii];
	}
	//indexing

	T at(int i, int j, int k)
	{
		uint ti = min(max(i, 0), (int)_nx - 1), tj = min(max(j, 0), (int)_ny - 1), tk = min(max(k, 0), (int)_nz - 1);
		//if(i>=0 && i<_nx && j>=0 && j<_ny && k>=0 && k<_nz )
		{
			uint I = ti >> 3, J = tj >> 3, K = tk >> 3;
			uint ii = ti & 7, jj = tj & 7, kk = tk & 7;
			uint idx = ((K*_blockx*_blocky + J*_blockx + I) << 9) + (kk << 6) + (jj << 3) + ii;
			return (_data->getPtr()[idx]);
			//return (index[K][J][I])[(((kk<<3)+jj)<<3)+ii];
		}
		//else
		//{
		//	return (T)0;
		//}
	}
	inline float lerp(float a, float b, float c)
	{
		return (1.0 - c)*(float)a + c*(float)b;
	}
	//sampling
	T sample_linear(double world_x, double world_y, double world_z)
	{
		double grid_x = world_x / _hx + _ox;
		double grid_y = world_y / _hy + _oy;
		double grid_z = world_z / _hz + _oz;

		int grid_i = (int)floor(grid_x);
		int grid_j = (int)floor(grid_y);
		int grid_k = (int)floor(grid_z);

		double cx = grid_x - (double)grid_i;
		double cy = grid_y - (double)grid_j;
		double cz = grid_z - (double)grid_k;

		float v1 = lerp(lerp((float)at(grid_i, grid_j, grid_k), (float)at(grid_i + 1, grid_j, grid_k), cx),
			lerp((float)at(grid_i, grid_j + 1, grid_k), (float)at(grid_i + 1, grid_j + 1, grid_k), cx),
			cy);
		float v2 = lerp(lerp((float)at(grid_i, grid_j, grid_k + 1), (float)at(grid_i + 1, grid_j, grid_k + 1), cx),
			lerp((float)at(grid_i, grid_j + 1, grid_k + 1), (float)at(grid_i + 1, grid_j + 1, grid_k + 1), cx),
			cy);
		return (T)(lerp(v1, v2, cz));

	}

	T sample_cubic(double world_x, double world_y, double world_z)
	{
		double gx = world_x / _hx + _ox;
		double gy = world_y / _hy + _oy;
		double gz = world_z / _hz + _oz;

		int i = (int)floor(gx);
		int j = (int)floor(gy);
		int k = (int)floor(gz);

		double x = gx - (double)i;
		double y = gy - (double)j;
		double z = gz - (double)k;

		float fn = cubic_interp((float)(at(i - 1, j - 1, k - 1)), (float)(at(i, j - 1, k - 1)),
			(float)(at(i + 1, j - 1, k - 1)), (float)(at(i + 2, j - 1, k - 1)), x);
		float f0 = cubic_interp((float)(at(i - 1, j, k - 1)), (float)(at(i, j, k - 1)),
			(float)(at(i + 1, j, k - 1)), (float)(at(i + 2, j, k - 1)), x);
		float f1 = cubic_interp((float)(at(i - 1, j + 1, k - 1)), (float)(at(i, j + 1, k - 1)),
			(float)(at(i + 1, j + 1, k - 1)), (float)(at(i + 2, j + 1, k - 1)), x);
		float f2 = cubic_interp((float)(at(i - 1, j + 2, k - 1)), (float)(at(i, j + 2, k - 1)),
			(float)(at(i + 1, j + 2, k - 1)), (float)(at(i + 2, j + 2, k - 1)), x);

		float fzn = cubic_interp(fn, f0, f1, f2, y);

		fn = cubic_interp((float)(at(i - 1, j - 1, k)), (float)(at(i, j - 1, k)),
			(float)(at(i + 1, j - 1, k)), (float)(at(i + 2, j - 1, k)), x);
		f0 = cubic_interp((float)(at(i - 1, j, k)), (float)(at(i, j, k)),
			(float)(at(i + 1, j, k)), (float)(at(i + 2, j, k)), x);
		f1 = cubic_interp((float)(at(i - 1, j + 1, k)), (float)(at(i, j + 1, k)),
			(float)(at(i + 1, j + 1, k)), (float)(at(i + 2, j + 1, k)), x);
		f2 = cubic_interp((float)(at(i - 1, j + 2, k)), (float)(at(i, j + 2, k)),
			(float)(at(i + 1, j + 2, k)), (float)(at(i + 2, j + 2, k)), x);

		float fz0 = cubic_interp(fn, f0, f1, f2, y);


		fn = cubic_interp((float)(at(i - 1, j - 1, k + 1)), (float)(at(i, j - 1, k + 1)),
			(float)(at(i + 1, j - 1, k + 1)), (float)(at(i + 2, j - 1, k + 1)), x);
		f0 = cubic_interp((float)(at(i - 1, j, k + 1)), (float)(at(i, j, k + 1)),
			(float)(at(i + 1, j, k + 1)), (float)(at(i + 2, j, k + 1)), x);
		f1 = cubic_interp((float)(at(i - 1, j + 1, k + 1)), (float)(at(i, j + 1, k + 1)),
			(float)(at(i + 1, j + 1, k + 1)), (float)(at(i + 2, j + 1, k + 1)), x);
		f2 = cubic_interp((float)(at(i - 1, j + 2, k + 1)), (float)(at(i, j + 2, k + 1)),
			(float)(at(i + 1, j + 2, k + 1)), (float)(at(i + 2, j + 2, k + 1)), x);

		float fz1 = cubic_interp(fn, f0, f1, f2, y);


		fn = cubic_interp((float)(at(i - 1, j - 1, k + 2)), (float)(at(i, j - 1, k + 2)),
			(float)(at(i + 1, j - 1, k + 2)), (float)(at(i + 2, j - 1, k + 2)), x);
		f0 = cubic_interp((float)(at(i - 1, j, k + 2)), (float)(at(i, j, k + 2)),
			(float)(at(i + 1, j, k + 2)), (float)(at(i + 2, j, k + 2)), x);
		f1 = cubic_interp((float)(at(i - 1, j + 1, k + 2)), (float)(at(i, j + 1, k + 2)),
			(float)(at(i + 1, j + 1, k + 2)), (float)(at(i + 2, j + 1, k + 2)), x);
		f2 = cubic_interp((float)(at(i - 1, j + 2, k + 2)), (float)(at(i, j + 2, k + 2)),
			(float)(at(i + 1, j + 2, k + 2)), (float)(at(i + 2, j + 2, k + 2)), x);

		float fz2 = cubic_interp(fn, f0, f1, f2, y);


		float res = cubic_interp(fzn, fz0, fz1, fz2, z);


		/*float v0,v1,v2,v3,v4,v5,v6,v7;
		sample_cube(world_x,world_y,world_z,v0,v1,v2,v3,v4,v5,v6,v7);
		float minv = min(min(min(min(min(min(min(v0,v1),v2),v3),v4),v5),v6),v7);
		float maxv = max(max(max(max(max(max(max(v0,v1),v2),v3),v4),v5),v6),v7);

		return max(min(maxv,res),minv);*/
		return res;

	}


	T sample_cube_lerp(double world_x, double world_y, double world_z, float &v0, float &v1, float &v2, float &v3, float &v4, float &v5, float &v6, float &v7)
	{
		double grid_x = world_x / _hx + _ox;
		double grid_y = world_y / _hy + _oy;
		double grid_z = world_z / _hz + _oz;

		int grid_i = (int)floor(grid_x);
		int grid_j = (int)floor(grid_y);
		int grid_k = (int)floor(grid_z);

		double cx = grid_x - (double)grid_i;
		double cy = grid_y - (double)grid_j;
		double cz = grid_z - (double)grid_k;

		v0 = (float)at(grid_i, grid_j, grid_k); v1 = (float)at(grid_i + 1, grid_j, grid_k);
		v2 = (float)at(grid_i, grid_j + 1, grid_k); v3 = (float)at(grid_i + 1, grid_j + 1, grid_k);
		v4 = (float)at(grid_i, grid_j, grid_k + 1); v5 = (float)at(grid_i + 1, grid_j, grid_k + 1);
		v6 = (float)at(grid_i, grid_j + 1, grid_k + 1); v7 = (float)at(grid_i + 1, grid_j + 1, grid_k + 1);

		float iv1 = lerp(lerp((float)at(grid_i, grid_j, grid_k), (float)at(grid_i + 1, grid_j, grid_k), cx),
			lerp((float)at(grid_i, grid_j + 1, grid_k), (float)at(grid_i + 1, grid_j + 1, grid_k), cx),
			cy);
		float iv2 = lerp(lerp((float)at(grid_i, grid_j, grid_k + 1), (float)at(grid_i + 1, grid_j, grid_k + 1), cx),
			lerp((float)at(grid_i, grid_j + 1, grid_k + 1), (float)at(grid_i + 1, grid_j + 1, grid_k + 1), cx),
			cy);
		return (T)(lerp(iv1, iv2, cz));

	}
	void sample_cube(double world_x, double world_y, double world_z, float &v0, float &v1, float &v2, float &v3, float &v4, float &v5, float &v6, float &v7)
	{
		double grid_x = world_x / _hx + _ox;
		double grid_y = world_y / _hy + _oy;
		double grid_z = world_z / _hz + _oz;

		int grid_i = (int)floor(grid_x);
		int grid_j = (int)floor(grid_y);
		int grid_k = (int)floor(grid_z);

		double cx = grid_x - (double)grid_i;
		double cy = grid_y - (double)grid_j;
		double cz = grid_z - (double)grid_k;

		v0 = (float)at(grid_i, grid_j, grid_k); v1 = (float)at(grid_i + 1, grid_j, grid_k);
		v2 = (float)at(grid_i, grid_j + 1, grid_k); v3 = (float)at(grid_i + 1, grid_j + 1, grid_k);
		v4 = (float)at(grid_i, grid_j, grid_k + 1); v5 = (float)at(grid_i + 1, grid_j, grid_k + 1);
		v6 = (float)at(grid_i, grid_j + 1, grid_k + 1); v7 = (float)at(grid_i + 1, grid_j + 1, grid_k + 1);

	}
};

typedef Buffer3D<float> buffer3Df;
typedef Buffer3D<double> buffer3Dd;
typedef Buffer3D<char> buffer3Dc;
typedef Buffer3D<int> buffer3Di;
typedef Buffer3D<uint> buffer3Dui;

#endif