// +-------------------------------------------------------------------------
// | cork.h
// | 
// | Author: Gilbert Bernstein
// +-------------------------------------------------------------------------
// | COPYRIGHT:
// |    Copyright Gilbert Bernstein 2013
// |    See the included COPYRIGHT file for further details.
// |    
// |    This file is part of the Cork library.
// |
// |    Cork is free software: you can redistribute it and/or modify
// |    it under the terms of the GNU Lesser General Public License as
// |    published by the Free Software Foundation, either version 3 of
// |    the License, or (at your option) any later version.
// |
// |    Cork is distributed in the hope that it will be useful,
// |    but WITHOUT ANY WARRANTY; without even the implied warranty of
// |    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// |    GNU Lesser General Public License for more details.
// |
// |    You should have received a copy 
// |    of the GNU Lesser General Public License
// |    along with Cork.  If not, see <http://www.gnu.org/licenses/>.
// +-------------------------------------------------------------------------
#pragma once

#include <string>

// WIN DLL export
// #ifdef _Windows ?
#if defined(_WINDLL)
#define CORK_IMP_EXP __declspec(dllexport)
#elif defined(CORK_STATIC)
#define CORK_IMP_EXP
#else
#define CORK_IMP_EXP __declspec(dllimport)
#endif


namespace Cork
{
	#ifndef uint
	typedef unsigned int uint;
	#endif

	// if a mesh is taken as input, the client must manage the memory
	// if a mesh is given as output, please use the provided
	// function to free the allocated memory.
	struct CORK_IMP_EXP CorkTriMesh
	{
		uint    n_triangles;
		uint    n_vertices;
		uint    *triangles;
		float   *vertices;

		CorkTriMesh()
			: triangles(0)
			, vertices(0)
		{}
	};

	CORK_IMP_EXP void freeCorkTriMesh(CorkTriMesh *mesh);

	// the inputs to Boolean operations must be "solid":
	//  -   closed (aka. watertight; see comment at bottom)
	//  -   non-self-intersecting
	// additionally, inputs should use a counter-clockwise convention
	// for triangle facing.  If the triangles are presented in clockwise
	// orientation, the object is interpreted as its unbounded complement

	// This function will test whether or not a mesh is solid
	CORK_IMP_EXP bool isSolid(CorkTriMesh mesh);

	// Boolean operations follow
	// result = A U B
	CORK_IMP_EXP void computeUnion(CorkTriMesh in0, CorkTriMesh in1, CorkTriMesh *out);

	// result = A - B
	CORK_IMP_EXP void computeDifference(CorkTriMesh in0, CorkTriMesh in1, CorkTriMesh *out);

	// result = A ^ B
	CORK_IMP_EXP void computeIntersection(CorkTriMesh in0, CorkTriMesh in1, CorkTriMesh *out);

	// result = A XOR B
	CORK_IMP_EXP void computeSymmetricDifference(
							CorkTriMesh in0, CorkTriMesh in1, CorkTriMesh *out);

	// Not a Boolean operation, but related:
	//  No portion of either surface is deleted.  However, the
	//  curve of intersection between the two surfaces is made explicit,
	//  such that the two surfaces are now connected.
	CORK_IMP_EXP void resolveIntersections(CorkTriMesh in0, CorkTriMesh in1, CorkTriMesh *out);


	class CORK_IMP_EXP Problem
	{
	public:
		/// Structure to store error codes
		struct ErrCode
		{
			static const int SUCCEED            = 0;
			static const int FAILED             = 1;
			static const int BAD_ALLOC          = 2;
			static const int UNKNOWN_OP         = 3;
			static const int WRONG_INPUT_SIZE   = 4;
		};

		/// Structure to share data with Cork
		struct DataExchanger
		{
			double * m_vertices;// pointer to the array of double corresponding to vertices
			unsigned m_numVertices;// number of vertices
			unsigned m_vertexSize;// size of a vertices must be almost 3, only the first three value of each vertices are used (as x, y and z).
			unsigned * m_triangles;// pointer to the array of unsigned corresponding to triangle, (one triangle correspond to three unsigned)
			unsigned m_numTriangles;// number of triangles
		};

		/// Define available boolean operations.
		enum Operator
		{
			Union,
			Difference,
			Intersection
		};

		/// Default constructor.
		Problem();

		/// Destructor.
		~Problem();

		int compute(const DataExchanger & A, Operator P, const DataExchanger & B);

		/// Get the number of vertices in result mesh.
		unsigned getNumVertices() const { return m_resNumVertices; }

		/// Get the number of triangles in result mesh.
		unsigned getNumTriangles() const { return m_resNumTriangles; }

		/// Get string corresponding to error returned by compute() or fillResult().
		std::string getErrorStr() const { return m_errorStr; }

		/// Fill resulting internal mesh to an already allocated DataExcanger.
		int fillResult(DataExchanger & p_result);

	private:
		/// clear internal mesh and result
		void clear();

		unsigned m_resNumTriangles;// number of mesh in result, initialized to zero
		unsigned m_resNumVertices;// number of mesh in result, initialized to zero
		void * m_mesh;
		std::string m_errorStr;
	};// class Problem
}// namespace Cork