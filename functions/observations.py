""" Functions to create observation points and cross sections for D3D-FM """

import os


class ObservationPoint():
    def __init__(self, name: str="", x: float=0, y: float=0):
        self.name = name
        self.x = x
        self.y = y


class ObservationCrossSection():
    def __init__(self, name: str="", x: list=[0, 1], y: list=[0, 1]):
        assert isinstance(x, list)
        assert isinstance(y, list)
        assert len(x) > 1
        assert len(x) == len(y)

        self.name = name
        self.x = x
        self.y = y
        self.n = len(x)


def write_observations(data: list, filename: str=None):
    """ Write observation points and observation cross sections to files 
    Input:
        data:       List of observation points and observation cross sections
        filename:   Name of file to write (extensions are added automatically)
    """
    if filename is None:
        filename = os.path.dirname(os.path.realpath(__file__)) + "/observations"
    if filename.endswith(".xyn"):
        filename.replace(".xyn", "")
    if filename.endswith("_crs.pli"):
        filename.replace("_crs.pli", "")

    filename_points = filename + ".xyn"
    filename_sections = filename + "_crs.pli"
    
    assert isinstance(data, list)

    # Sort points and cross sections
    points = []
    sections = []

    for element in data:
        if isinstance(element, ObservationPoint):
            points.append(element)
        elif isinstance(element, ObservationCrossSection):
            sections.append(element)
        else:
            print(f"'{element}' is not an ObservationPoint or ObservationCrossSection")
    
    # Process points
    if len(points) > 0:
        with open(filename_points, "w") as file:
            for element in points:
                _name = element.name.replace("\'", "")
                file.write(f"{element.x} \t{element.y} \t'{_name}'\n")

    # Process cross sections
    if len(sections) > 0:
        with open(filename_sections, "w") as file:
            for element in sections:
                _name = element.name.replace("\'", "")
                file.write(f"{_name}\n")
                file.write(f"{element.n} \t2\n")
                for i in range(element.n):
                    file.write(f"\t{element.x[i]} \t{element.y[i]} \n")
                file.write(f"\n")
    
    return


if __name__ == "__main__":
    obs_0 = ObservationPoint(name="Center", x=0, y=0)
    obs_1 = ObservationPoint(name="Point", x=1, y=1)
    obs_2 = ObservationCrossSection(name="Line", x=[-2, 2], y=[0, 0])
    obs_3 = ObservationCrossSection(name="Curve", x=[-2, 0, 0], y=[0, 0, 2])
    obs_4 = ObservationCrossSection(name="Square", x=[-4, -4 ,4, 4], y=[-4, 4, 4, -4])

    write_observations([obs_0, obs_1, obs_2, obs_3, obs_4])
