"""Module for linear algebra operations without numpy."""

import math

def angle(vector1, vector2):
    """Angle (in radians) between vectors vector1 and vector2."""
    return math.acos(dot(vector1, vector2) /
                     (normvector(vector1) * normvector(vector2)))


def cross(vector1, vector2):
    """Cross product of vectors vector1 and vector2."""
    return [vector1[1]*vector2[2] - vector1[2]*vector2[1],
            vector1[2]*vector2[0] - vector1[0]*vector2[2],
            vector1[0]*vector2[1] - vector1[1]*vector2[0]]


def dot(vector1, vector2):
    """Dot product of vectors vector1 and vector2."""
    return sum((a*b) for a, b in zip(vector1, vector2))


def euclideandistance(point1, point2):
    """Euclidean distance between points point1 and point2."""
    return math.sqrt((point1[0] - point2[0])**2 +
                     (point1[1] - point2[1])**2 +
                     (point1[2] - point2[2])**2)


def normvector(vector):
    """Norm of vector."""
    return math.sqrt(dot(vector, vector))


def normalizevector(vector):
    """Normalize vector."""
    norm = normvector(vector)
    return [vector[0] / norm,
            vector[1] / norm,
            vector[2] / norm]


def project(vector1, vector2):
    """Project vector vector1 along vector vector2."""
    vector2norm = normalizevector(vector2)
    vector1component = dot(vector1, vector2norm)
    return [vector1component * vector2norm[0],
            vector1component * vector2norm[1],
            vector1component * vector2norm[2]]


def tangent(vector, normal):
    """Tangential component of vector for given normal."""
    return [vector[0] - dot(vector, normal) * normal[0],
            vector[1] - dot(vector, normal) * normal[1],
            vector[2] - dot(vector, normal) * normal[2]]
