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

class Vertex 
{
public:
	glm::vec3 position;
	glm::mat4 Q; // QEM at each vertex
	std::set<int> adjacentFaces; // Index of adjacent face
	bool deleteFlag = false;
};

class Face
{
public:
	std::vector<int> v;
	glm::vec3 n;
	float d;
	bool isBoundary = false;
};

class Edge
{
public:
	int v1, v2;
	double cost;
	glm::vec3 optimalPos; // Optimal vertex position with minimum cost
	int version = 0;

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
	void computeErrorMetrics(int vertexIndex);
	void computeEdgeCost(Edge& edge);
	void updateArrays();
	bool validEdge(const Edge& edge);
	bool isFlipped(const Edge& edge, const std::set<int>& commonFaces);
	void simplify();

public:
	std::vector<Vertex> vertices;
	std::vector<Face> faceVertices;
	std::vector<int> directedEdge;
	std::vector<int> oppositeHalf;

	std::vector<glm::vec3> normals;
	std::vector<glm::vec2> texCoords;

	std::priority_queue<Edge> edgeHeap;
	std::map<std::pair<int, int>, int> versionControl;

	std::vector<int> indexMap;

	// BoundingBox
	glm::vec3 minBBX;
	glm::vec3 maxBBX;

	char fileName[64] = "";

	const char* SimplifyMode[2] = { "Optimal Position", "Middle Position"};
	int currentMode = 0;
	float simplifyRatio = 0.01f;
	float t = 0.2f;

	bool preserveTexture = false;
	bool preserveBoundary = false;
	bool preventInversion = false;
	bool enableAggregation = false;
	bool enableTexture = false;
	bool hasNormal = false;
	bool hasTexture = false;

	//float error = 0.0f;
};