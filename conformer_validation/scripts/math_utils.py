import math
import numpy as np

def bat_value(coords):
    if coords.shape[0] == 2:
        v = MathUtils.Distance(*coords)
    elif coords.shape[0] == 3:
        v = MathUtils.Angle(*coords, units='radians')
    elif coords.shape[0] == 4:
        v = MathUtils.Dihedral(*coords, units='radians')
    else:
        raise Exception(
            'The length of coords is wrong. Coords: {}'.format(coords))
    return v

class MathUtils(object):
    @staticmethod
    def Vector(p1, p2):
        """calculate distance

        Args:
            p1 (numpy.array): point 1
            p2 (numpy.array): point 2

        Returns:
            numpy.array: vector V_12
        """
        return p2 - p1

    @staticmethod
    def Distance(p1, p2):
        """calculate distance

        Args:
            p1 (numpy.array): point 1
            p2 (numpy.array): point 2

        Returns:
            float: distance_12
        """
        return MathUtils.distance(p1, p2)

    @staticmethod
    def Angle(p1, p2, p3, units="degrees"):
        """calculate angle

        Args:
            p1 (numpy.array): point 1
            p2 (numpy.array): point 2
            p3 (numpy.array): point 3

        Returns:
            float: angle_123
        """
        return MathUtils.angle(p1-p2, p3-p2, units=units)

    @staticmethod
    def Dihedral(p1, p2, p3, p4, units="degrees"):
        """calculate dihedral

        Args:
            p1 (numpy.array): point 1
            p2 (numpy.array): point 2
            p3 (numpy.array): point 3
            p4 (numpy.array): point 4

        Returns:
            float: dihedral_1234
        """
        return MathUtils.dihedral(p2-p1, p3-p2, p4-p3, units=units)


    @staticmethod
    def distance(v1, v2):
        """calculate distance

        Args:
            v1 (numpy.array): point 1
            v2 (numpy.array): point 2

        Returns:
            float: distance
        """
        return np.linalg.norm(v1 - v2)


    @staticmethod
    def angle(v1, v2, units="degrees"):
        """Calculates the angle between two vectors.

        copied from pymatgen

        Args:
            v1: Vector 1
            v2: Vector 2
            units: "degrees" or "radians". Defaults to "degrees".

        Returns:
            Angle between them in degrees.
        """
        d = np.dot(v1, v2) / np.linalg.norm(v1) / np.linalg.norm(v2)
        d = min(d, 1)
        d = max(d, -1)
        angle = math.acos(d)
        if units == "degrees":
            return math.degrees(angle)
        elif units == "radians":
            return angle
        else:
            raise ValueError("Invalid units {}".format(units)) # pragma: no cover

    @staticmethod
    def dihedral(v1, v2, v3, units='degrees'):
        """Calculates the dihedral between three vectors.

        Args:
            v1: Vector 1
            v2: Vector 2
            v3: Vector 3
            units: "degrees" or "radians". Defaults to "degrees".

        Returns:
            Angle between them in degrees.
        """
        v12 = np.cross(v1, v2)
        v23 = np.cross(v2, v3)
        dihe = math.atan2(np.dot(v12, v3) * np.linalg.norm(v2), np.dot(v12, v23))
        if units == "degrees":
            return math.degrees(dihe)
        elif units == "radians":
            return dihe
        else:
            raise ValueError("Invalid units {}".format(units)) # pragma: no cover

def __calcrms(y1_coord, y2_coord):
    y1c = y1_coord - sum(y1_coord) / len(y1_coord)
    y2c = y2_coord - sum(y2_coord) / len(y2_coord)
    return rmsdxyz().kabsch_rmsd(y1c, y2c)


def calcrms(mol1, mol2, remove_h=True):
    import numpy as np

    coords1, coords2 = mol1.GetConformer(0).GetPositions() , mol2.GetConformer(0).GetPositions() 

    if remove_h:
        elements = []
        for a1, a2 in zip(mol1.GetAtoms(), mol2.GetAtoms()):
            assert a1.GetAtomicNum() == a2.GetAtomicNum()
            elements.append(a1.GetAtomicNum())

    if not remove_h:
        return __calcrms(coords1, coords2)

    c1, c2 = [], []
    for i, ele in enumerate(elements):
        if ele == 1:
            continue
        c1.append(coords1[i])
        c2.append(coords2[i])

    assert c1 and c2
    return __calcrms(np.array(c1), np.array(c2))

class rmsdxyz(object):
    def __init__(self):
        pass

    def kabsch_rmsd(self,P, Q):
        """
        Rotate matrix P unto Q and calculate the RMSD
        """
        P = self.rotate(P, Q)
        return self.rmsd(P, Q)


    def rotate(self,P, Q):
        """
        Rotate matrix P unto matrix Q using Kabsch algorithm
        """
        U = self.kabsch(P, Q)

        # Rotate P
        P = np.dot(P, U)
        return P


    def kabsch(self,P, Q):
        """
        The optimal rotation matrix U is calculated and then used to rotate matrix
        P unto matrix Q so the minimum root-mean-square deviation (RMSD) can be
        calculated.

        Using the Kabsch algorithm with two sets of paired point P and Q,
        centered around the center-of-mass.
        Each vector set is represented as an NxD matrix, where D is the
        the dimension of the space.

        The algorithm works in three steps:
        - a translation of P and Q
        - the computation of a covariance matrix C
        - computation of the optimal rotation matrix U

        http://en.wikipedia.org/wiki/Kabsch_algorithm

        Parameters:
        P -- (N, number of points)x(D, dimension) matrix
        Q -- (N, number of points)x(D, dimension) matrix

        Returns:
        U -- Rotation matrix

        """

        # Computation of the covariance matrix
        C = np.dot(np.transpose(P), Q)

        # Computation of the optimal rotation matrix
        # This can be done using singular value decomposition (SVD)
        # Getting the sign of the det(V)*(W) to decide
        # whether we need to correct our rotation matrix to ensure a
        # right-handed coordinate system.
        # And finally calculating the optimal rotation matrix U
        # see http://en.wikipedia.org/wiki/Kabsch_algorithm
        V, S, W = np.linalg.svd(C)
        d = (np.linalg.det(V) * np.linalg.det(W)) < 0.0

        if d:
            S[-1] = -S[-1]
            V[:, -1] = -V[:, -1]

        # Create Rotation matrix U
        U = np.dot(V, W)

        return U


    def centroid(self,X):
        """
        Calculate the centroid from a vectorset X
        """
        C = sum(X)/len(X)
        return C


    def rmsd(self,V, W):
        """
        Calculate Root-mean-square deviation from two sets of vectors V and W.
        """
        D = len(V[0])
        N = len(V)
        rmsd = 0.0
        for v, w in zip(V, W):
            rmsd += sum([(v[i]-w[i])**2.0 for i in range(D)])
        return np.sqrt(rmsd/N)


    def write_coordinates(self,atoms, V):
        """
        Print coordinates V
        """
        N, D = V.shape

        print(str(N))
        print()

        for i in range(N):
            line = "{0:2s} {1:15.8f} {2:15.8f} {3:15.8f}".format(atoms[i], V[i, 0], V[i, 1], V[i, 2])
            print(line)

    def get_coordinates(self,filename, fmt, ignore_hydrogens):
        """
        Get coordinates from filename.

        """
        if fmt == "xyz":
            return self.get_coordinates_xyz(filename, ignore_hydrogens)
        if fmt == "pdb":
            return self.get_coordinates_pdb(filename, ignore_hydrogens)
        exit("Could not recognize file format: {:s}".format(fmt))


    def get_coordinates_pdb(self,filename, ignore_hydrogens):
        """
        Get coordinates from the first chain in a pdb file
        and return a vectorset with all the coordinates.

        """
        # PDB files tend to be a bit of a mess. The x, y and z coordinates
        # are supposed to be in column 31-38, 39-46 and 47-54, but this is not always the case.
        # Because of this the three first columns containing a decimal is used.
        # Since the format doesn't require a space between columns, we use the above
        # column indices as a fallback.
        x_column = None
        V = []
        # Same with atoms and atom naming. The most robust way to do this is probably
        # to assume that the atomtype is given in column 3.
        atoms = []
        with open(filename) as f:
            lines = f.readlines()
            for line in lines:
                if line.startswith("TER") or line.startswith("END"):
                    break
                if line.startswith("ATOM"):
                    tokens = line.split()
                    # Try to get the atomtype
                    try:
                        atom = tokens[2][0]
                        if ignore_hydrogens and atom == "H":
                            continue
                        elif atom in ["H", "C", "N", "O", "S", "P"]:
                            atoms.append(atom)
                    except:
                            exit("Error parsing atomtype for the following line: \n%s" % line)

                    if x_column == None:
                        try:
                            # look for x column
                            for i, x in enumerate(tokens):
                                if "." in x and "." in tokens[i+1] and "." in tokens[i+2]:
                                    x_column = i
                                    break
                        except IndexError:
                            exit("Error parsing coordinates for the following line: \n%s" % line)
                    # Try to read the coordinates
                    try:
                        V.append(np.asarray(tokens[x_column:x_column+3],dtype=float))
                    except:
                        # If that doesn't work, use hardcoded indices
                        try:
                            x = line[30:38]
                            y = line[38:46]
                            z = line[46:54]
                            V.append(np.asarray([x,y,z],dtype=float))
                        except:
                            exit("Error parsing input for the following line: \n%s" % line)


        V = np.asarray(V)
        return atoms, V


    def get_coordinates_xyz(self,filename, ignore_hydrogens):
        """
        Get coordinates from a filename.xyz and return a vectorset with all the
        coordinates.

        This function has been written to parse XYZ files, but can easily be
        written to parse others.

        """
        import re
        
        f = open(filename, 'r')
        V = []
        atoms = []
        n_atoms = 0
        lines_read = 0

        # Read the first line to obtain the number of atoms to read
        try:
            n_atoms = int(next(f))
        except ValueError:
            exit("Could not obtain the number of atoms in the .xyz file.")

        # Skip the title line
        next(f)

        # Use the number of atoms to not read beyond the end of a file
        for line in f:

            if lines_read == n_atoms:
                break

            atom = re.findall(r'[a-zA-Z]+', line)[0]
            numbers = re.findall(r'[-]?\d+\.\d*(?:[Ee][-\+]\d+)?', line)
            numbers = [float(number) for number in numbers]

            # ignore hydrogens
            if ignore_hydrogens and atom.lower() == "h":
                continue

            # The numbers are not valid unless we obtain exacly three
            if len(numbers) == 3:
                V.append(np.array(numbers))
                atoms.append(atom)
            else:
                exit("Reading the .xyz file failed in line {0}. Please check the format.".format(lines_read + 2))

            lines_read += 1

        f.close()
        V = np.array(V)
        return atoms, V
