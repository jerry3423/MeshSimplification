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
			Vertex vertex;
			s >> vertex.position[0] >> vertex.position[1] >> vertex.position[2];
			vertices.push_back(vertex);
		}
		else if (pre == "vn")
		{
			// Load normals
			hasNormal = true;
			glm::vec3 normal;
			s >> normal[0] >> normal[1] >> normal[2];
			normals.push_back(normal);
		}
		else if (pre == "vt")
		{
			// Load texture coordinates
			hasTexture = true;
			glm::vec2 texCoord;
			s >> texCoord[0] >> texCoord[1];
			texCoords.push_back(texCoord);
		}
		else if (pre == "f")
		{
			Face face;
			face.v.resize(3);

			// Load indices of position/ texture coordinate/ normal
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
	// Add "-simp" to indicate it's a simplified mesh
	std::string filePath = baseFilePath + filename + "-simp" + extension;
	ofstream outfile(filePath);

	for (const auto& vertex : vertices)
	{
		outfile << "v " << fixed << setprecision(6) << vertex.position[0] << " " << vertex.position[1] << " " << vertex.position[2] << "\n";
	}
	// Write normals to the file if they exist
	if (hasNormal)
	{
		for (const auto& normal : normals)
		{
			outfile << "vn " << fixed << setprecision(6) << normal[0] << " " << normal[1] << " " << normal[2] << "\n";
		}
	}
	// Write texture coordinates to the file if they exist
	if (hasTexture)
	{
		for (const auto& texCoord : texCoords)
		{
			outfile << "vt " << fixed << setprecision(6) << texCoord[0] << " " << texCoord[1] << "\n";
		}
	}
	// Write face vertex indices to the file
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
			// Construct a half edge from vertex vFrom to vertex vTo
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
		// Create the opposite half edge
		std::pair<int, int> oppositeHalfEdge(halfEdge.second, halfEdge.first);
		// If opposite half edge exists, map them to each other
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
	// Build the directed edge structure for the mesh
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

			// Get the index of the next half edge
			hEdgeIndex = oppositeHalf[hEdgeIndex];
			if (hEdgeIndex != -1)
			{
				hEdgeIndex = 3 * (hEdgeIndex / 3) + (hEdgeIndex % 3 + 1) % 3;
			}
		} while (hEdgeIndex != -1 && hEdgeIndex != i);
	}

	// Mark boundary vertices
	for (size_t i = 0; i < directedEdge.size(); i++)
	{
		if (oppositeHalf[i] == -1)
		{
			vertices[directedEdge[i]].isBoundary = true;
		}
	}

	// Compute the plane equation parameters for each face
	for (size_t i = 0; i < faceVertices.size(); i++)
	{
		computeFaceParameter(i);
	}

	// Compute Q for each vertex
	for (size_t i = 0; i < vertices.size(); i++)
	{
		computeQuadricError(i);
	}

	// Build the priority queue for edges based on their collapse cost
	while (!edgeHeap.empty())
	{
		edgeHeap.pop();
	}
	std::map<std::pair<int, int>, bool> processedEdges;
	for (size_t i = 0; i < directedEdge.size(); i++)
	{
		// Get indices of two vertices
		int v1 = directedEdge[i];
		int v2 = directedEdge[3 * (i / 3) + (i % 3 + 1) % 3];

		// Avoid edge has been added to edge queue
		if (v1 > v2) std::swap(v1, v2);

		auto edgePair = std::make_pair(v1, v2);
		// Continue if this edge has been added
		if (processedEdges.find(edgePair)!=processedEdges.end()) continue;
		processedEdges[edgePair] = true;

		// Create edge and push to priority queue
		Edge edge;
		edge.v1 = v1;
		edge.v2 = v2;
		edge.version = 0;
		computeEdgeCost(edge);
		edgeHeap.push(edge);
		versionControl[edgePair] = 0; // Default edge version to 0
	}

	// Additional processing for close vertices aggregation if enabled
	if (enableAggregation)
	{
		// Find all vertices whose distance from the current vertex is less than t
		std::vector<std::vector<int>> closeVertices(vertices.size());
		for (size_t i = 0; i < vertices.size(); i++) {
			for (size_t j = i + 1; j < vertices.size(); j++) {
				if ((vertices[i].position - vertices[j].position).norm() < t) {
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
				versionControl[edgePair] = 0; // Default edge version to 0
			}
		}
	}

}

void Mesh::computeFaceParameter(int faceIndex)
{
	// Get three vertices
	Eigen::Vector3f v1 = vertices[faceVertices[faceIndex].v[0]].position;
	Eigen::Vector3f v2 = vertices[faceVertices[faceIndex].v[1]].position;
	Eigen::Vector3f v3 = vertices[faceVertices[faceIndex].v[2]].position;

	// Calculate the face normal
	Eigen::Vector3f faceNormal = (v2 - v1).cross(v3 - v1).normalized();

	// Update the face normal and distance d in the face structure
	faceVertices[faceIndex].n = faceNormal;
	faceVertices[faceIndex].d = -faceNormal.dot(v1);
}

void Mesh::computeQuadricError(int vertexIndex)
{
	// Initialize matrix Q
	vertices[vertexIndex].Q.setZero();

	// Iterate over all faces adjacent to this vertex
	for (auto it = vertices[vertexIndex].adjacentFaces.begin(); it != vertices[vertexIndex].adjacentFaces.end(); it++)
	{
		Face& face = faceVertices[*it];

		// Calculate the fundamental error matrix components A, b, c for the face
		Eigen::Matrix3f A = face.n * face.n.transpose();
		Eigen::Vector3f b = face.d * face.n;
		float c = face.d * face.d;

		// Construct the 4x4 error matrix for this face
		Eigen::Matrix4f Kp;
		Kp << A, b,
			b.transpose(), c;

		// Add the error matrix for this face to the vertex's quadric error matrix
		vertices[vertexIndex].Q += Kp;
	}
}

void Mesh::computeEdgeCost(Edge& edge)
{
	Eigen::Matrix4f Q = vertices[edge.v1].Q + vertices[edge.v2].Q;

	// Extract A, b, and c from matrix Q
	Eigen::Matrix3f A = Q.topLeftCorner<3, 3>();
	Eigen::Vector3f b = Q.topRightCorner<3, 1>();
	float c = Q(3, 3);

	// Attempt to invert matrix A
	Eigen::Matrix3f Ainv;
	bool isAInvertible;
	A.computeInverseWithCheck(Ainv, isAInvertible, 1e-2);
	// Optimal position mode
	if (currentMode == 0)
	{
		// If A is invertible, calculate the optimal position and cost
		if (isAInvertible)
		{
			edge.optimalPos = -Ainv * b;
			edge.cost = b.dot(edge.optimalPos) + c;
		}
		else
		{
			// If A is not invertible, use the midpoint as the optimal position
			edge.optimalPos = (vertices[edge.v1].position + vertices[edge.v2].position) * 0.5f;
			Eigen::Vector4f v(edge.optimalPos(0), edge.optimalPos(1), edge.optimalPos(2), 1.0f);
			edge.cost = v.transpose() * Q * v;
		}
	}
	else if (currentMode == 1)
	{
		// Use the midpoint as the optimal position
		edge.optimalPos = (vertices[edge.v1].position + vertices[edge.v2].position) * 0.5f;
		Eigen::Vector4f v(edge.optimalPos(0), edge.optimalPos(1), edge.optimalPos(2), 1.0f);
		edge.cost = v.transpose() * Q * v;
	}
}

void Mesh::updateArrays()
{
	// Mark vertices with no adjacent faces for deletion
	for (size_t i = 0; i < vertices.size(); i++)
	{
		if (vertices[i].adjacentFaces.size() == 0)
		{
			vertices[i].deleteFlag = true;
		}
	}
	// Prepare new arrays for updated mesh data
	std::vector<Vertex> newV;
	std::vector<Face> newF;
	std::vector<glm::vec2> newTex;
	std::vector<glm::vec3> newNorm;

	// Map to keep track of new indices of vertices
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
			indexMap[i] = -1; // Mark deleted vertices with -1
		}
	}

	for (size_t i = 0; i < faceVertices.size(); i++) {
		// Check if all vertices of the face are still present
		if (!vertices[faceVertices[i].v[0]].deleteFlag && !vertices[faceVertices[i].v[1]].deleteFlag && !vertices[faceVertices[i].v[2]].deleteFlag) {
			Face updatedFace = faceVertices[i];
			for (int j = 0; j < 3; j++) {
				updatedFace.v[j] = indexMap[updatedFace.v[j]];
			}
			newF.push_back(updatedFace);
		}
	}

	// Replace old arrays with the new ones
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
	// Perserve boundaries
	if ((vertices[v1].isBoundary || vertices[v2].isBoundary) && preserveBoundary)
	{
		return false;
	}
	return true;
}

bool Mesh::checkFaceInversion(const Edge& edge, const std::set<int>& commonFaces)
{
	int v1 = edge.v1;
	int v2 = edge.v2;
	// Get all adjacent faces
	std::set<int> facesToCheck;
	facesToCheck.insert(vertices[v1].adjacentFaces.begin(), vertices[v1].adjacentFaces.end());
	facesToCheck.insert(vertices[v2].adjacentFaces.begin(), vertices[v2].adjacentFaces.end());
	// Erase common adjacent faces
	for (int faceIndex : commonFaces)
	{
		facesToCheck.erase(faceIndex);
	}

	for (int faceIndex : facesToCheck)
	{
		// Calculate old normal
		Face face = faceVertices[faceIndex];
		Eigen::Vector3f oldNormal = (vertices[face.v[1]].position - vertices[face.v[0]].position).cross(vertices[face.v[2]].position - vertices[face.v[0]].position).normalized();

		// Get three vertex positions if edge collapse happen
		Eigen::Vector3f newPositions[3];
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
		Eigen::Vector3f newNormal = (newPositions[1] - newPositions[0]).cross(newPositions[2] - newPositions[0]).normalized();

		float d = std::cos(40 * 3.1415926 / 180.0);
		// Check if the direction of the normal has changed significantly
		if (oldNormal.dot(newNormal) < d)
		{
			return true;
		}
	}

	return false;

}

void Mesh::simplify()
{
	// Calculate the target number of faces based on the simplify ratio
	int face_num = faceVertices.size();
	size_t count = 0;
	while (count < (1 - simplifyRatio) * face_num)
	{
		// Select the edge with the minimum cost that is valid for collapsing
		while (!edgeHeap.empty() && !validEdge(edgeHeap.top()))
		{
			edgeHeap.pop();
		}
		if (edgeHeap.empty())
		{
			break;
		}
		// Get the edge with the minimum cost
		Edge edge = edgeHeap.top();
		edgeHeap.pop();

		int v1 = edge.v1;
		int v2 = edge.v2;

		// Get v1 and v2 common adjacent faces
		std::set<int> commonAdjacentFaces;
		for (int faceIndex : vertices[v1].adjacentFaces)
		{
			if (vertices[v2].adjacentFaces.find(faceIndex) != vertices[v2].adjacentFaces.end())
			{
				commonAdjacentFaces.insert(faceIndex);
			}
		}
		// Skip edges that would create non-manifold geometry
		if (commonAdjacentFaces.size() > 2)
		{
			continue;
		}

		// Prevent face inversion
		if (preventInversion)
		{
			if (checkFaceInversion(edge, commonAdjacentFaces))
			{
				continue;
			}
		}

		// Update v1 position to optimal opsition
		vertices[v1].position = edge.optimalPos;
		vertices[v1].Q = vertices[v1].Q + vertices[v2].Q;

		// Remove collapsed faces that are adjacent to both v1 and v2
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

		// Mark v2 as deleted
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

		// Rebuild the edges connected to v1 and add them to the priority queue
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
		// Increment the count based on how many faces were collapsed
		count += commonAdjacentFaces.size();
	}

	updateArrays();

}
