#pragma once

#include<iostream>
#include<vector>
#include<string>
#include<map>
#include<set>
#include<queue>
#include<functional>
#include<unordered_map>

#include<glm/glm.hpp>
#include<Eigen/Dense>
#include<Eigen/LU>

class Vertex 
{
public:
	Eigen::Vector3f position; // vertex position
	Eigen::Matrix4f Q; // QEM at each vertex
	std::set<int> adjacentFaces; // Index of adjacent face
	bool deleteFlag = false; // Whether the vertex is deleted during edge collapse operation
	bool isBoundary = false; // Whether the vertex is boundary
};

class Face
{
public:
	std::vector<int> v; // Indices of vertices that form this face
	Eigen::Vector3f n; // Normal vector of face
	float d;
};

class Edge
{
public:
	int v1, v2; // Two indices of vertices on an edge
	double cost; // Cost of doing edge collapse operation
	Eigen::Vector3f optimalPos; // Optimal vertex position with minimum cost
	int version = 0; // version of the edge, used in priority queue for update check

	bool operator<(const Edge& other) const {
		return cost > other.cost;
	}
};


class Mesh {

public:
	Mesh();
	void readFile(const char* filename);
	void writeFile(const char* filename);
	void clear();
	void buildDirectedEdgeStructure();
	void preProcessData();
	void computeFaceParameter(int faceIndex);
	void computeQuadricError(int vertexIndex);
	void computeEdgeCost(Edge& edge);
	void updateArrays();
	bool validEdge(const Edge& edge);
	bool checkFaceInversion(const Edge& edge, const std::set<int>& commonFaces);
	void simplify();

public:
	std::vector<Vertex> vertices;
	std::vector<Face> faceVertices;
	std::vector<int> directedEdge;
	std::vector<int> oppositeHalf;

	std::vector<glm::vec3> normals;
	std::vector<glm::vec2> texCoords;

	std::priority_queue<Edge> edgeHeap; // Priority queue for edges based on cost
	std::map<std::pair<int, int>, int> versionControl;

	std::vector<int> indexMap; // Map from old to new vertex indices

	// BoundingBox
	glm::vec3 minBBX;
	glm::vec3 maxBBX;

	char fileName[64] = ""; // Model file name

	const char* SimplifyMode[2] = { "Optimal Position", "Middle Position"};
	int currentMode = 0;
	float simplifyRatio = 0.01f;
	float t = 0.2f; //Distance threshold

	bool preserveTexture = false;
	bool preserveBoundary = true;
	bool preventInversion = true;
	bool enableAggregation = false;
	bool enableTexture = false;
	bool hasNormal = false;
	bool hasTexture = false;

};