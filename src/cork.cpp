// +-------------------------------------------------------------------------
// | cork.cpp
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
#include "cork.h"

#include "mesh/mesh.h"
#include "mesh/corkMesh.h"

namespace Cork
{
	void freeCorkTriMesh(CorkTriMesh *mesh)
	{
		if (mesh->triangles)
			delete[] mesh->triangles;
		if (mesh->vertices)
			delete[] mesh->vertices;
		mesh->n_triangles = 0;
		mesh->n_vertices = 0;
	}

	void corkTriMesh2CorkMesh(
		CorkTriMesh in,
		CorkMesh *mesh_out
	) {
		RawCorkMesh raw;
		raw.vertices.resize(in.n_vertices);
		raw.triangles.resize(in.n_triangles);
		if(in.n_vertices == 0 || in.n_triangles == 0) {
			CORK_ERROR("empty mesh input to Cork routine.");
			*mesh_out = CorkMesh(raw);
			return;
		}
    
		uint max_ref_idx = 0;
		for(uint i=0; i<in.n_triangles; i++) {
			raw.triangles[i].a = in.triangles[3*i+0];
			raw.triangles[i].b = in.triangles[3*i+1];
			raw.triangles[i].c = in.triangles[3*i+2];
			max_ref_idx = std::max(
							std::max(max_ref_idx,
									 in.triangles[3*i+0]),
							std::max(in.triangles[3*i+1],
									 in.triangles[3*i+2])
						  );
		}
		if(max_ref_idx > in.n_vertices) {
			CORK_ERROR("mesh input to Cork routine has an out of range reference "
				  "to a vertex.");
			raw.vertices.clear();
			raw.triangles.clear();
			*mesh_out = CorkMesh(raw);
			return;
		}
    
		for(uint i=0; i<in.n_vertices; i++) {
			raw.vertices[i].pos.x = in.vertices[3*i+0];
			raw.vertices[i].pos.y = in.vertices[3*i+1];
			raw.vertices[i].pos.z = in.vertices[3*i+2];
		}
    
		*mesh_out = CorkMesh(raw);
	}
	void corkMesh2CorkTriMesh(
		CorkMesh *mesh_in,
		CorkTriMesh *out
	) {
		RawCorkMesh raw = mesh_in->raw();
    
		out->n_triangles = raw.triangles.size();
		out->n_vertices  = raw.vertices.size();
    
		out->triangles = new uint[(out->n_triangles) * 3];
		out->vertices  = new float[(out->n_vertices) * 3];
    
		for(uint i=0; i<out->n_triangles; i++) {
			(out->triangles)[3*i+0] = raw.triangles[i].a;
			(out->triangles)[3*i+1] = raw.triangles[i].b;
			(out->triangles)[3*i+2] = raw.triangles[i].c;
		}
    
		for(uint i=0; i<out->n_vertices; i++)
		{
			(out->vertices)[3*i+0] = static_cast<float>(raw.vertices[i].pos.x);
			(out->vertices)[3*i+1] = static_cast<float>(raw.vertices[i].pos.y);
			(out->vertices)[3*i+2] = static_cast<float>(raw.vertices[i].pos.z);
		}
	}


	bool isSolid(CorkTriMesh cmesh)
	{
		CorkMesh mesh;
		corkTriMesh2CorkMesh(cmesh, &mesh);
    
		bool solid = true;
    
		if(mesh.isSelfIntersecting()) {
			CORK_ERROR("isSolid() was given a self-intersecting mesh");
			solid = false;
		}
    
		if(!mesh.isClosed()) {
			CORK_ERROR("isSolid() was given a non-closed mesh");
			solid = false;
		}
    
		return solid;
	}

	void computeUnion(
		CorkTriMesh in0, CorkTriMesh in1, CorkTriMesh *out
	) {
		CorkMesh cmIn0, cmIn1;
		corkTriMesh2CorkMesh(in0, &cmIn0);
		corkTriMesh2CorkMesh(in1, &cmIn1);
    
		cmIn0.boolUnion(cmIn1);
    
		corkMesh2CorkTriMesh(&cmIn0, out);
	}

	void computeDifference(
		CorkTriMesh in0, CorkTriMesh in1, CorkTriMesh *out
	) {
		CorkMesh cmIn0, cmIn1;
		corkTriMesh2CorkMesh(in0, &cmIn0);
		corkTriMesh2CorkMesh(in1, &cmIn1);
    
		cmIn0.boolDiff(cmIn1);
    
		corkMesh2CorkTriMesh(&cmIn0, out);
	}

	void computeIntersection(
		CorkTriMesh in0, CorkTriMesh in1, CorkTriMesh *out
	) {
		CorkMesh cmIn0, cmIn1;
		corkTriMesh2CorkMesh(in0, &cmIn0);
		corkTriMesh2CorkMesh(in1, &cmIn1);
    
		cmIn0.boolIsct(cmIn1);
    
		corkMesh2CorkTriMesh(&cmIn0, out);
	}

	void computeSymmetricDifference(
		CorkTriMesh in0, CorkTriMesh in1, CorkTriMesh *out
	) {
		CorkMesh cmIn0, cmIn1;
		corkTriMesh2CorkMesh(in0, &cmIn0);
		corkTriMesh2CorkMesh(in1, &cmIn1);
    
		cmIn0.boolXor(cmIn1);
    
		corkMesh2CorkTriMesh(&cmIn0, out);
	}

	void resolveIntersections(
		CorkTriMesh in0, CorkTriMesh in1, CorkTriMesh *out
	) {
		CorkMesh cmIn0, cmIn1;
		corkTriMesh2CorkMesh(in0, &cmIn0);
		corkTriMesh2CorkMesh(in1, &cmIn1);
    
		cmIn0.disjointUnion(cmIn1);
		cmIn0.resolveIntersections();
    
		corkMesh2CorkTriMesh(&cmIn0, out);
	}

	///////////////////////////////////////////////////////////////////////////
	//
	//              class problem implementation
	//
	///////////////////////////////////////////////////////////////////////////

	//-----------------------------------------------------------------------
	Problem::Problem() :
	//-----------------------------------------------------------------------
		m_mesh(NULL),
		m_resNumTriangles(0),
		m_resNumVertices(0),
		m_errorStr("no error, or error not set")
	{
	}

	//-----------------------------------------------------------------------
	Problem::~Problem()
	//-----------------------------------------------------------------------
	{
		clear();
	}

	/// Function to fill a CorkMesh from a DataExchanger
	bool convert(CorkMesh * p_mesh, const Problem::DataExchanger & A)
	{
		std::vector<CorkVertex>& outVerts = p_mesh->vertices();
		std::vector<CorkMesh::Tri>& outTriangles = p_mesh->triangles();

		try
		{
			outTriangles.resize(A.m_numTriangles);
			outVerts.resize(A.m_numVertices);
		}
		catch(std::bad_alloc)
		{
			//not enough memory
			return false;
		}

		// fill out the vertices
		double * ptr = A.m_vertices;
		for (size_t i=0; i < A.m_numVertices; i++)
		{
			outVerts[i].pos.x = static_cast<double>( *(ptr++) );
			outVerts[i].pos.y = static_cast<double>( *(ptr++) );
			outVerts[i].pos.z = static_cast<double>( *(ptr++) );
			// ptr += A.m_vertexSize - 3;
			for (unsigned i = 3; i < A.m_vertexSize; ++i)
				++ptr;
			
		}

		// fill out the triangles
		for (size_t i=0; i < A.m_numTriangles; i++)
		{
			//DGM: strangely, Cork seems to duplicate the (a,b,c) information...
			outTriangles[i].data.a  = outTriangles[i].a = static_cast<int>(A.m_triangles[ 3 * i    ]);
			outTriangles[i].data.b  = outTriangles[i].b = static_cast<int>(A.m_triangles[ 3 * i + 1]);
			outTriangles[i].data.c  = outTriangles[i].c = static_cast<int>(A.m_triangles[ 3 * i + 2]);
		}

		return true;
	}// convert(...)

	//-----------------------------------------------------------------------
	void Problem::clear()
	//-----------------------------------------------------------------------
	{
		m_resNumTriangles = m_resNumVertices = 0;
		if (m_mesh)
		{
			delete m_mesh;
			m_mesh = NULL;
		}
	}

	//-----------------------------------------------------------------------
	int Problem::compute(const DataExchanger & A,
						 Operator P,
						 const DataExchanger & B)
	//-----------------------------------------------------------------------
	{
		// clear previous result
		clear();

		if (A.m_numVertices == 0 || A.m_numTriangles == 0 || A.m_vertexSize < 3 ||
			B.m_numVertices == 0 || B.m_numTriangles == 0 || B.m_vertexSize < 3)
		{
			//empty input mesh or wrong usage
			m_errorStr = "[apply] wrong input size";
			return ErrCode::WRONG_INPUT_SIZE;
		}

		m_mesh = new CorkMesh();

		CorkMesh * corkMeshA = static_cast<CorkMesh *>(m_mesh);
		CorkMesh corkMeshB;

		if (convert(corkMeshA, A) && convert(&corkMeshB, B))
		{
			switch (P)
			{
			case Intersection:
				corkMeshA->boolIsct(corkMeshB);
				break;
			case Union:
				corkMeshA->boolUnion(corkMeshB);
				break;
			case Difference:
				corkMeshA->boolDiff(corkMeshB);
				break;
			default:
				m_errorStr = "[apply] unknown operator";
				clear();
				return ErrCode::UNKNOWN_OP;
			}

			m_resNumTriangles = corkMeshA->numTris();
			m_resNumVertices = corkMeshA->numVerts();

			return true;
		}
		else
		{
			m_errorStr = "[apply::convert] bad allocation";
			clear();
			return ErrCode::BAD_ALLOC;
		}

		return ErrCode::SUCCEED;
	}// compute(...)

	//-----------------------------------------------------------------------
	int Problem::fillResult(DataExchanger & p_result)
	//-----------------------------------------------------------------------
	{
		CorkMesh * corkMesh = static_cast<CorkMesh *>(m_mesh);
		if ( corkMesh->numTris() != p_result.m_numTriangles ||
			 corkMesh->numTris() != p_result.m_numTriangles )
		{
			m_errorStr = "[fillResult] wrong input data size";
			return ErrCode::WRONG_INPUT_SIZE;
		}

		const std::vector<CorkVertex>& inVerts = corkMesh->vertices();
		const std::vector<CorkMesh::Tri>& inTriangles = corkMesh->triangles();

		// fill out the vertices
		double * outVerts = p_result.m_vertices;
		for (size_t i=0; i<inVerts.size(); i++)
		{
			*(outVerts++) = inVerts[i].pos.x;
			*(outVerts++) = inVerts[i].pos.y;
			*(outVerts++) = inVerts[i].pos.z;

			// outVerts += A.p_result.m_vertexSize - 3;
			for (unsigned j = 3; j < p_result.m_vertexSize; ++j)
				++outVerts;
		}

		// fill out the triangles
		unsigned * outTriangles = p_result.m_triangles;
		for (size_t i=0; i<inTriangles.size(); i++)
		{
			*(outTriangles++) = inTriangles[i].a;
			*(outTriangles++) = inTriangles[i].b;
			*(outTriangles++) = inTriangles[i].c;
		}

		return ErrCode::SUCCEED;
	}// fillResult(...)
}// namespace Cork