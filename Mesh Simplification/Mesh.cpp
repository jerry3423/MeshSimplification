#include"Mesh.h"

#include<iostream>
#include<fstream>
#include<sstream>
#include<iomanip>
#include<map>
#include<algorithm>
#include<iterator>

using namespace std;

std::string baseFilePath = "../models/";
std::string extension = ".obj";

// Init mesh
Mesh::Mesh()
{

}

void Mesh::readFile(const char* filename)
{
	std::string filePath = baseFilePath + filename + extension;
	ifstream infile(filePath);

	string line;
	while (getline(infile, line)) {
		istringstream s(line);
		string pre;
		s >> pre;

		size_t nIndex = 0, tIndex = 0;
		if (pre == "v")
		{
			// Load vertex position
			// v 0.732802 -0.238102 1.137424
			Vertex vertex;
			s >> vertex.position[0] >> vertex.position[1] >> vertex.position[2];
			vertices.push_back(vertex);
		}
		else if (pre == "vn")
		{
			// Load normals
			// vn -0.873412 0.234567 0.594623
			hasNormal = true;
			glm::vec3 normal;
			s >> normal[0] >> normal[1] >> normal[2];
			normals.push_back(normal);
		}
		else if (pre == "vt")
		{
			// Load texture coordinates
			// vt 0.900000 1.000000
			hasTexture = true;
			glm::vec2 texCoord;
			s >> texCoord[0] >> texCoord[1];
			texCoords.push_back(texCoord);
		}
		else if (pre == "f")
		{
			Face face;
			face.v.resize(3);

			// Load vertex attributes                     (v/vn/vt)
			// If texture coordinates and normals exists: f 1/1/1 2/2/2 3/3/3
			// If only texture coordinates exists:        f 1//1 2//2 3//3
			for (int i = 0; i < 3; ++i) {
				std::string vertexAttri;
				s >> vertexAttri;
				std::istringstream vertexStream(vertexAttri);
				std::string v, vt, vn;
				std::getline(vertexStream, v, '/');
				std::getline(vertexStream, vt, '/');
				std::getline(vertexStream, vn);
				face.v[i] = std::stoi(v) - 1;
			}
			faceVertices.push_back(face);
		}
	}

	infile.close();
}

void Mesh::writeFile(const char* filename)
{
	std::string filePath = baseFilePath + filename + "-simp" + extension;
	ofstream outfile(filePath);

	for (const auto& vertex : vertices)
	{
		outfile << "v " << fixed << setprecision(6) << vertex.position[0] << " " << vertex.position[1] << " " << vertex.position[2] << "\n";
	}
	if (hasNormal)
	{
		for (const auto& normal : normals)
		{
			outfile << "vn " << fixed << setprecision(6) << normal[0] << " " << normal[1] << " " << normal[2] << "\n";
		}
	}
	if (hasTexture)
	{
		for (const auto& texCoord : texCoords)
		{
			outfile << "vt " << fixed << setprecision(6) << texCoord[0] << " " << texCoord[1] << "\n";
		}
	}
	for (const auto& face : faceVertices)
	{
		outfile << "f " << fixed << setprecision(6) << face.v[0] + 1 << " " << face.v[1] + 1 << " " << face.v[2] + 1 << "\n";
	}

	outfile.close();
}

void Mesh::clear()
{
	vertices.resize(0);
	faceVertices.resize(0);
	directedEdge.resize(0);
	oppositeHalf.resize(0);
	texCoords.resize(0);
	normals.resize(0);
	hasTexture = false;
	hasNormal = false;
}

void Mesh::buildDirectedEdgeStructure()
{
	directedEdge.clear();
	oppositeHalf.clear();
	std::map<std::pair<int, int>, int> halfEdgeMap;

	// Iterate over all the triangles
	for (size_t i = 0; i < faceVertices.size(); i++)
	{
		// Iterate over the three half edge of the triangle
		for (int j = 0; j < 3; j++)
		{
			// Build a half edge from v0 to v1
			int vFrom = faceVertices[i].v[j];
			int vTo = faceVertices[i].v[(j + 1) % 3];
			std::pair<int, int> halfEdge(vFrom, vTo);

			// Build first directed edge
			if (halfEdgeMap.find(halfEdge) == halfEdgeMap.end())
			{
				int halfEdgeIndex = directedEdge.size();
				directedEdge.push_back(vFrom);
				halfEdgeMap[halfEdge] = halfEdgeIndex;
			}
		}
	}

	// Build opposite half edge
	oppositeHalf.resize(directedEdge.size(), -1);
	for (const auto& edgePair : halfEdgeMap)
	{
		std::pair<int, int> halfEdge = edgePair.first;
		std::pair<int, int> oppositeHalfEdge(halfEdge.second, halfEdge.first);
		if (halfEdgeMap.find(oppositeHalfEdge) != halfEdgeMap.end())
		{
			int halfEdgeIndex = halfEdgeMap[halfEdge];
			int oppositeIndex = halfEdgeMap[oppositeHalfEdge];
			oppositeHalf[halfEdgeIndex] = oppositeIndex;
		}
	}

}

void Mesh::preProcessData()
{
	buildDirectedEdgeStructure();

	for (int i = 0; i < vertices.size(); i++)
	{
		vertices[i].adjacentFaces.clear();
	}
	std::vector<bool> visitedHalfEdges(directedEdge.size(), false);

	// Get vertex's adjacent faces
	for (size_t i = 0; i < directedEdge.size(); i++)
	{
		if (visitedHalfEdges[i]) continue;

		int vFromIndex = directedEdge[i];
		int hEdgeIndex = i;

		do {
			visitedHalfEdges[hEdgeIndex] = true;
			vertices[vFromIndex].adjacentFaces.insert(hEdgeIndex / 3);

			// Get next halfedge Index
			hEdgeIndex = oppositeHalf[hEdgeIndex];
			if (hEdgeIndex != -1)
			{
				hEdgeIndex = 3 * (hEdgeIndex / 3) + (hEdgeIndex % 3 + 1) % 3;
			}
		} while (hEdgeIndex != -1 && hEdgeIndex != i);
	}

	// Get boundary vertex
	for (size_t i = 0; i < directedEdge.size(); i++)
	{
		if (oppositeHalf[i] == -1)
		{
			vertices[directedEdge[i]].isBoundary = true;
		}
	}

	// Compute parameters of each face equation: ax + by + cz + d = 0
	for (size_t i = 0; i < faceVertices.size(); i++)
	{
		computeFaceParameter(i);
	}

	// Compute Q for each vertex
	for (size_t i = 0; i < vertices.size(); i++)
	{
		computeErrorMetrics(i);
	}
	// Build edge queue
	// Clear queue first
	while (!edgeHeap.empty())
	{
		edgeHeap.pop();
	}
	std::map<std::pair<int, int>, bool> processedEdges;
	for (size_t i = 0; i < directedEdge.size(); i++)
	{
		int v1 = directedEdge[i];
		int v2 = directedEdge[3 * (i / 3) + (i % 3 + 1) % 3];

		if (v1 > v2) std::swap(v1, v2);

		auto edgePair = std::make_pair(v1, v2);
		// Continue if this edge has been added
		if (processedEdges.find(edgePair)!=processedEdges.end()) continue;
		processedEdges[edgePair] = true;

		//  Create edge and push to priority queue
		Edge edge;
		edge.v1 = v1;
		edge.v2 = v2;
		edge.version = 0;
		computeEdgeCost(edge);
		edgeHeap.push(edge);
		versionControl[edgePair] = 0;
	}

	if (enableAggregation)
	{
		// Find all vertices whose distance from the current vertex is less than t
		std::vector<std::vector<int>> closeVertices(vertices.size());
		for (size_t i = 0; i < vertices.size(); i++) {
			for (size_t j = i + 1; j < vertices.size(); j++) {
				if (glm::distance(vertices[i].position, vertices[j].position) < t) {
					closeVertices[i].push_back(j);
					closeVertices[j].push_back(i);
				}
			}
		}

		// Build these vertex pair and push to edge queue
		for (size_t i = 0; i < closeVertices.size(); i++)
		{
			for (int j : closeVertices[i])
			{
				int v1 = i;
				int v2 = j;

				if (v1 > v2) std::swap(v1, v2);
				// Continue if this edge has been added
				auto edgePair = std::make_pair(v1, v2);
				if (processedEdges.find(edgePair) != processedEdges.end()) continue;
				processedEdges[edgePair] = true;

				Edge edge;
				edge.v1 = v1;
				edge.v2 = v2;
				edge.version = 0;
				computeEdgeCost(edge);
				edgeHeap.push(edge);
				versionControl[edgePair] = 0;
			}
		}
	}

}

void Mesh::computeFaceParameter(int faceIndex)
{
	//Get three vertices
	glm::vec3 v1 = vertices[faceVertices[faceIndex].v[0]].position;
	glm::vec3 v2 = vertices[faceVertices[faceIndex].v[1]].position;
	glm::vec3 v3 = vertices[faceVertices[faceIndex].v[2]].position;

	glm::vec3 faceNormal = glm::cross(v2 - v1, v3 - v1);

	faceVertices[faceIndex].a = faceNormal.x;
	faceVertices[faceIndex].b = faceNormal.y;
	faceVertices[faceIndex].c = faceNormal.z;
	faceVertices[faceIndex].d = -glm::dot(faceNormal, v1);
}

void Mesh::computeErrorMetrics(int vertexIndex)
{
	// Initial matrix Q
	vertices[vertexIndex].Q = glm::mat4(0.0f);

	// Iterate certex's adjacent faces and compute Q
	for (auto it = vertices[vertexIndex].adjacentFaces.begin(); it != vertices[vertexIndex].adjacentFaces.end(); it++)
	{
		Face face = faceVertices[*it];
		glm::mat4 Kp(
			face.a * face.a, face.a * face.b, face.a * face.c, face.a * face.d,
			face.a * face.b, face.b * face.b, face.b * face.c, face.b * face.d,
			face.a * face.c, face.c * face.b, face.c * face.c, face.c * face.d,
			face.d * face.a, face.d * face.b, face.d * face.c, face.d * face.d
		);

		vertices[vertexIndex].Q += Kp;
	}
	// Evaluate geometry error
	glm::vec4 v(vertices[vertexIndex].position, 1.0f);
	d_squared += glm::dot(v, vertices[vertexIndex].Q * v);

	// If boundary preservation
	if (preserveBoundary)
	{
		if (vertices[vertexIndex].isBoundary)
		{
			vertices[vertexIndex].Q *= 1000;
		}
	}
}

void Mesh::computeEdgeCost(Edge& edge)
{
	glm::mat4 QP = vertices[edge.v1].Q + vertices[edge.v2].Q;
	glm::mat4 Q = vertices[edge.v1].Q + vertices[edge.v2].Q;
	Q[3][0] = 0.0f;
	Q[3][1] = 0.0f;
	Q[3][2] = 0.0f;
	Q[3][3] = 1.0f;

	glm::vec4 p;
	if (std::abs(glm::determinant(Q)) > 1e-6)
	{
		glm::mat4 b = glm::inverse(Q);
		p = b * glm::vec4(0.0f, 0.0f, 0.0f, 1.0f);
	}
	else {
		p = glm::vec4((vertices[edge.v1].position + vertices[edge.v2].position) * 0.5f, 1.0f);
	}

	edge.optimalPos = glm::vec3(p.x, p.y, p.z);
	edge.cost = glm::dot(p, QP * p);
}

void Mesh::updateArrays()
{
	// Clean vertices data
	for (size_t i = 0; i < vertices.size(); i++)
	{
		if (vertices[i].adjacentFaces.size() == 0)
		{
			vertices[i].deleteFlag = true;
		}
	}
	// Update origin vertex array and face array
	std::vector<Vertex> newV;
	std::vector<Face> newF;
	std::vector<glm::vec2> newTex;
	std::vector<glm::vec3> newNorm;


	indexMap.resize(vertices.size());
	int simplifiedindex = 0;
	for (size_t i = 0; i < vertices.size(); i++)
	{
		// If vertex has not been deleted, update index and push to a new vertex array
		if (!vertices[i].deleteFlag)
		{
			indexMap[i] = simplifiedindex;
			simplifiedindex++;
			newV.push_back(vertices[i]);
			if (hasTexture)
			{
				newTex.push_back(texCoords[i]);
			}
			if (hasNormal)
			{
				newNorm.push_back(normals[i]);
			}
		}
		else {
			indexMap[i] = -1;
		}
	}

	for (size_t i = 0; i < faceVertices.size(); i++) {
		// If three vertices in a face have not been deleted, push correct vertex index to a new face array
		if (!vertices[faceVertices[i].v[0]].deleteFlag && !vertices[faceVertices[i].v[1]].deleteFlag && !vertices[faceVertices[i].v[2]].deleteFlag) {
			Face updatedFace = faceVertices[i];
			for (int j = 0; j < 3; j++) {
				updatedFace.v[j] = indexMap[updatedFace.v[j]];
			}
			newF.push_back(updatedFace);
		}
	}

	vertices = std::move(newV);
	faceVertices = std::move(newF);
	texCoords = std::move(newTex);
	normals = std::move(newNorm);
}

bool Mesh::validEdge(const Edge& edge)
{
	int v1 = edge.v1;
	int v2 = edge.v2;

	// If v1 or v2 has been deleted, return false
	if (vertices[v1].deleteFlag || vertices[v2].deleteFlag)
	{
		return false;
	}

	// Check if edge v1 - v2 is the latest
	auto it = versionControl.find(std::make_pair(edge.v1, edge.v2));
	if (it==versionControl.end() || edge.version != it->second)
	{
		return false;
	}

	if (edge.cost > 1e4)
	{
		return false;
	}

	return true;
}

bool Mesh::isFlipped(const Edge& edge, const std::set<int>& commonFaces)
{
	int v1 = edge.v1;
	int v2 = edge.v2;
	// Get all adjacent faces
	std::set<int> facesToCheck;
	facesToCheck.insert(vertices[v1].adjacentFaces.begin(), vertices[v1].adjacentFaces.end());
	facesToCheck.insert(vertices[v2].adjacentFaces.begin(), vertices[v2].adjacentFaces.end());

	for (int faceIndex : commonFaces)
	{
		facesToCheck.erase(faceIndex);
	}

	for (int faceIndex : facesToCheck)
	{
		// Calculate old normal
		Face face = faceVertices[faceIndex];
		glm::vec3 oldNormal = glm::normalize(glm::cross(vertices[face.v[1]].position - vertices[face.v[0]].position, 
			vertices[face.v[2]].position - vertices[face.v[0]].position));

		// Get three vertex postion if edge collapse happen
		glm::vec3 newPositions[3];
		for (int i = 0; i < 3; i++)
		{
			if (face.v[i] == v1 || face.v[i] == v2)
			{
				newPositions[i] = edge.optimalPos;
			}
			else
			{
				newPositions[i] = vertices[face.v[i]].position;
			}
		}

		// Calculate new normal
		glm::vec3 newNormal = glm::normalize(glm::cross(newPositions[1] - newPositions[0], newPositions[2] - newPositions[0]));

		if (glm::dot(oldNormal, newNormal) < 0.0f)
		{
			return true;
		}
	}

	return false;

}

bool Mesh::commonVerticesNum(const Edge& edge, int commonFacesSize)
{
	// Get v1's adjacent vertices
	std::set<int> v1adjacentVertices;
	for (int faceIndex : vertices[edge.v1].adjacentFaces)
	{
		for (int i = 0; i < 3; i++)
		{
			if (faceVertices[faceIndex].v[i] == edge.v1)
			{
				continue;
			}
			v1adjacentVertices.insert(faceVertices[faceIndex].v[i]);
		}
	}

	// Get v2's adjacent vertices
	std::set<int> v2adjacentVertices;
	for (int faceIndex : vertices[edge.v2].adjacentFaces)
	{
		for (int i = 0; i < 3; i++)
		{
			if (faceVertices[faceIndex].v[i] == edge.v2)
			{
				continue;
			}
			v2adjacentVertices.insert(faceVertices[faceIndex].v[i]);
		}
	}

	int count = 0;
	for (auto v1index : v1adjacentVertices)
	{
		if (v2adjacentVertices.find(v1index) != v2adjacentVertices.end())
		{
			count++;
		}
	}

	if (count != commonFacesSize)
	{
		return false;
	}

	return true;
}

void Mesh::simplify()
{
	// Compute Q for all initial vertices

	// Loop through all edge pairs and build a priority edge queue

		// Compute optimal vertex position for each edge pair using Q1, Q2

		// Compute minimum cost for each edge pair

	// Repeat:
	
		// Select and collapse an edge pair(v1,v2) with minimum cost in the top of queue
		
		// Push new edge to queue
	//-----------------------------------------------------------------------------------

	int face_num = faceVertices.size();
	size_t count = 0;
	// Select and collapse edge, v1 - v2. 
	// Move v1 and v2 to optimal position, replace v2 to v1, delete v2
	while (count < (1 - simplifyRatio) * face_num)
	{
		// Get a vaild edge with minumum cost
		while (!edgeHeap.empty() && !validEdge(edgeHeap.top()))
		{
			edgeHeap.pop();
		}
		if (edgeHeap.empty())
		{
			break;
		}
		Edge edge = edgeHeap.top();
		edgeHeap.pop();

		int v1 = edge.v1;
		int v2 = edge.v2;

		// Get v1 and v2 common adjacent face idnex
		std::set<int> commonAdjacentFaces;
		for (int faceIndex : vertices[v1].adjacentFaces)
		{
			if (vertices[v2].adjacentFaces.find(faceIndex) != vertices[v2].adjacentFaces.end())
			{
				commonAdjacentFaces.insert(faceIndex);
			}
		}

		if (preventInversion)
		{
			if (isFlipped(edge, commonAdjacentFaces))
			{
				continue;
			}
		}

		/*if (!commonVerticesNum(edge, commonAdjacentFaces.size()))
		{
			continue;
		}*/
		// Update v1 position to optimal opsition
		vertices[v1].position = edge.optimalPos;
		vertices[v1].Q = vertices[v1].Q + vertices[v2].Q;

		// Erase vertex's adjacent faces that have edge v1 - v2
		for (int faceIndex : commonAdjacentFaces)
		{
			for (int i = 0; i < 3; i++)
			{
				vertices[faceVertices[faceIndex].v[i]].adjacentFaces.erase(faceIndex);
			}
		}

		// Update vertexIndex in each face from v2 to v1
		for (int faceIndex : vertices[v2].adjacentFaces)
		{
			for (int i = 0; i < 3; i++)
			{
				if (faceVertices[faceIndex].v[i] == v2)
				{
					faceVertices[faceIndex].v[i] = v1;
				}
			}
			// Merge v2's adjacent faces into v1
			vertices[v1].adjacentFaces.insert(faceIndex);
		}

		vertices[v2].deleteFlag = true;

		// Get v1's adjacent vertices
		std::set<int> adjacentVertices;
		for (int faceIndex : vertices[v1].adjacentFaces)
		{
			for (int i = 0; i < 3; i++)
			{
				if (faceVertices[faceIndex].v[i] == v1)
				{
					continue;
				}
				adjacentVertices.insert(faceVertices[faceIndex].v[i]);
			}
		}

		// Build new edge and push to edge queue
		for (int vertexIndex : adjacentVertices)
		{
			int vIndex = vertexIndex;
			Edge newedge;
			newedge.v1 = std::min(vIndex, v1);
			newedge.v2 = std::max(vIndex, v1);
			newedge.version = ++versionControl[std::make_pair(newedge.v1, newedge.v2)];
			computeEdgeCost(newedge);
			edgeHeap.push(newedge);
		}
		if (commonAdjacentFaces.size() == 1)
		{
			count++;
		}
		else if (commonAdjacentFaces.size() == 2)
		{
			count += 2;
		}
	}

	updateArrays();

}
