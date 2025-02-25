import os
import glob
import xml.etree.ElementTree as ET
from xml.dom import minidom
import sys
import subprocess

def open_file_in_default_viewer(file_path):
    if sys.platform.startswith("win"):
        os.startfile(file_path)
    elif sys.platform == "darwin":
        subprocess.run(["open", file_path])
    else:
        subprocess.run(["xdg-open", file_path])

def read_cap_info(project_dir):
    """
    Attempt to read 'cap_info.txt' from the given project directory.

    Returns:
      (areas, face_names): A tuple of two lists, sorted alphabetically by face_name.
         - areas: list of float (areas)
         - face_names: list of str (face names)

    Raises:
      FileNotFoundError if 'cap_info.txt' is not found.
      ValueError if a line is malformed (optional handling).
    """
    cap_file_path = os.path.join(project_dir, "cap_info.txt")

    # Raise an error if the file does not exist
    if not os.path.isfile(cap_file_path):
        raise FileNotFoundError(f"Could not find 'cap_info.txt' in project directory: {project_dir}")

    data = []  # Will store (face_name, area) tuples

    with open(cap_file_path, "r") as f:
        for line in f:
            line = line.strip()
            # Skip empty lines or comment lines (that start with '#')
            if not line or line.startswith('#'):
                continue

            parts = line.split()
            if len(parts) != 2:
                # You might choose to raise an error, skip, or handle differently
                # For now, let's raise an error
                raise ValueError(f"Malformed line in cap_info.txt: '{line}'")

            face_name, area_str = parts
            area = float(area_str)  # Convert area to float

            data.append((face_name, area))

    # Sort by face_name alphabetically
    data.sort(key=lambda x: x[0])

    # Separate into parallel lists
    face_names = [item[0] for item in data]
    areas = [item[1] for item in data]
    
    return areas, face_names

def prompt_BC(areas, face_names):
    """
    Prompts the user for:
      - total pressure drop,
      - left-right ratio,
      - inflow velocity for each face that does NOT start with 'cap_lpa' or 'cap_rpa'.
    
    Returns:
      A dictionary with keys:
        - "total_pressure_drop": float
        - "left_right_ratio": float
        - "inflows": dict mapping {face_name: inflow_velocity} for non-lpa/rpa faces
    """

    print("\n=== BC Prompt ===")
    
    # 1. Prompt for total pressure drop
    while True:
        try:
            total_pressure_drop = float(input("Enter the total pressure drop (mmHg): "))
            break
        except ValueError:
            print("Invalid entry. Please enter a numeric value for pressure drop.")

    # 2. Prompt for left-right ratio
    while True:
        try:
            left_ratio = float(input("Enter the pulmonary flow ratio (enter percentage flow to the left pulmonary artery): "))
            break
        except ValueError:
            print("Invalid entry. Please enter a numeric value for the ratio.")

    # 3. For each face NOT starting with cap_lpa or cap_rpa, prompt for inflow velocity
    inflows = []
    inflow_names = []
    for face, area in zip(face_names, areas):
        # Check if face name starts with either 'cap_lpa' or 'cap_rpa'
        if not (face.startswith("cap_lpa") or face.startswith("cap_rpa")):
            while True:
                try:
                    prompt_text = (f"Enter inflow velocity for face '{face}' "
                                   f"(area = {area}): ")
                    velocity = float(input(prompt_text))
                    inflows.append(velocity)
                    inflow_names.append(face)
                    break
                except ValueError:
                    print("Invalid entry. Please enter a numeric value for velocity.")
    return total_pressure_drop, left_ratio, inflows, inflow_names

def calculate_resistance(total_pressure_drop, inflows, left_ratio, areas, face_names):
    total_inflow = 0
    for v in inflows:
        total_inflow += v
        
    lpa = total_inflow * left_ratio/100
    rpa = total_inflow - lpa

    pressure = total_pressure_drop * 1333.2
    LR = pressure / lpa
    RR = pressure / rpa
    
    left_total_area = 0
    right_total_area = 0
    for a, face in zip(areas, face_names):
        if (face.startswith("cap_lpa")):
            left_total_area += a
        elif (face.startswith("cap_rpa")):
            right_total_area += a
    
    area_ratios = []
    cap_faces = []
    for a, face in zip(areas, face_names):
        if (face.startswith("cap_lpa")):
            area_ratios.append(a/left_total_area)
            cap_faces.append(face)
        elif (face.startswith("cap_rpa")):
            area_ratios.append(a/right_total_area)
            cap_faces.append(face)
    
    resistances = []
    for a, face in zip(area_ratios, cap_faces):
        if (face.startswith("cap_lpa")):
            resistances.append(pressure / (a * lpa))
        elif (face.startswith("cap_rpa")):
            resistances.append(pressure / (a * rpa))
    print(area_ratios)
    return resistances, cap_faces

    
    
def create_mesh_xml(project_dir, output_xml_path, resistances, cap_faces, inflows, inflow_names):
    """
    Create an XML file for your solver with:
      - <Add_face> + <Add_BC> only for .vtp files starting with 'cap'
        and also for 'walls_combined.vtp'.
      - Alphabetical order of these faces.
      - If face name starts with 'cap_lpa' or 'cap_rpa', use a Neu/Resistance BC.
      - Otherwise, use Dir/Steady/0.0 BC.
    """

    # 1. Check paths
    mesh_dir = os.path.join(project_dir, "mesh")
    mesh_surfaces_dir = os.path.join(mesh_dir, "mesh-surfaces")
    if not os.path.isdir(mesh_dir):
        raise FileNotFoundError(f"Mesh folder not found at: {mesh_dir}")
    if not os.path.isdir(mesh_surfaces_dir):
        raise FileNotFoundError(f"Mesh-surfaces folder not found at: {mesh_surfaces_dir}")

    # 2. Collect relevant .vtp files
    all_vtp_files = glob.glob(os.path.join(mesh_surfaces_dir, "*.vtp"))
    
    # Filter for those starting with 'cap' plus check for walls_combined.vtp
    filtered_vtp_files = []
    for vtp_file in all_vtp_files:
        filename = os.path.basename(vtp_file)  # e.g. 'cap_inlet.vtp'
        if filename.startswith("cap"):
            filtered_vtp_files.append(vtp_file)
    
    # Optionally add 'walls_combined.vtp' if present
    walls_combined_path = os.path.join(mesh_dir, "walls_combined.vtp")
    if os.path.isfile(walls_combined_path):
        filtered_vtp_files.append(walls_combined_path)

    # Sort alphabetically
    filtered_vtp_files.sort(key=lambda p: os.path.basename(p).lower())

    # 3. Root element
    root = ET.Element("svMultiPhysicsFile", {"version": "0.1"})
    root.tail = " "

    # 4. GeneralSimulationParameters
    gsp = ET.SubElement(root, "GeneralSimulationParameters")
    gsp.tail = " "

    ET.SubElement(gsp, "Continue_previous_simulation").text = "false"
    ET.SubElement(gsp, "Number_of_spatial_dimensions").text = "3"
    ET.SubElement(gsp, "Number_of_time_steps").text = "2"
    ET.SubElement(gsp, "Time_step_size").text = "0.005"
    ET.SubElement(gsp, "Spectral_radius_of_infinite_time_step").text = "0.50"
    ET.SubElement(gsp, "Searched_file_name_to_trigger_stop").text = "STOP_SIM"

    ET.SubElement(gsp, "Save_results_to_VTK_format").text = "1"
    ET.SubElement(gsp, "Name_prefix_of_saved_VTK_files").text = "result"
    ET.SubElement(gsp, "Increment_in_saving_VTK_files").text = "2"
    ET.SubElement(gsp, "Start_saving_after_time_step").text = "1"

    ET.SubElement(gsp, "Increment_in_saving_restart_files").text = "100"
    ET.SubElement(gsp, "Convert_BIN_to_VTK_format").text = "0"

    ET.SubElement(gsp, "Verbose").text = "1"
    ET.SubElement(gsp, "Warning").text = "0"
    ET.SubElement(gsp, "Debug").text = "0"

    # 5. <Add_mesh>
    add_mesh_elem = ET.SubElement(root, "Add_mesh", {"name": "msh"})
    add_mesh_elem.tail = " "

    mesh_file_path = ET.SubElement(add_mesh_elem, "Mesh_file_path")
    mesh_file_path.text = "mesh/mesh-complete.mesh.vtu"
    mesh_file_path.tail = " "
    
    # Create <Add_face> for each filtered file
    face_names = []
    for vtp_file in filtered_vtp_files:
        if (os.path.splitext(os.path.basename(vtp_file))[0] != "walls_combined"):
            face_name = os.path.splitext(os.path.basename(vtp_file))[0]
            add_face_elem = ET.SubElement(add_mesh_elem, "Add_face", {"name": face_name})
            face_file_path_elem = ET.SubElement(add_face_elem, "Face_file_path")
            face_file_path_elem.text = f"mesh/mesh-surfaces/{os.path.basename(vtp_file)}"
            add_face_elem.tail = " "
            face_names.append(face_name)
        else:
            face_name = os.path.splitext(os.path.basename(vtp_file))[0]
            add_face_elem = ET.SubElement(add_mesh_elem, "Add_face", {"name": face_name})
            face_file_path_elem = ET.SubElement(add_face_elem, "Face_file_path")
            face_file_path_elem.text = f"mesh/{os.path.basename(vtp_file)}"
            add_face_elem.tail = " "
            face_names.append(face_name)
    

    # 6. <Add_equation type="fluid">
    add_equation_elem = ET.SubElement(root, "Add_equation", {"type": "fluid"})

    ET.SubElement(add_equation_elem, "Coupled").text = "true"
    ET.SubElement(add_equation_elem, "Min_iterations").text = "3"
    ET.SubElement(add_equation_elem, "Max_iterations").text = "5"
    ET.SubElement(add_equation_elem, "Tolerance").text = "1e-11"
    ET.SubElement(add_equation_elem, "Backflow_stabilization_coefficient").text = "0.2"

    ET.SubElement(add_equation_elem, "Density").text = "1.06"
    viscosity_elem = ET.SubElement(add_equation_elem, "Viscosity", {"model": "Constant"})
    ET.SubElement(viscosity_elem, "Value").text = "0.04"

    # Output
    output_spatial_elem = ET.SubElement(add_equation_elem, "Output", {"type": "Spatial"})
    ET.SubElement(output_spatial_elem, "Velocity").text = "true"
    ET.SubElement(output_spatial_elem, "Pressure").text = "true"
    ET.SubElement(output_spatial_elem, "Traction").text = "true"
    ET.SubElement(output_spatial_elem, "Vorticity").text = "true"
    ET.SubElement(output_spatial_elem, "Divergence").text = "true"
    ET.SubElement(output_spatial_elem, "WSS").text = "true"

    output_bint_elem = ET.SubElement(add_equation_elem, "Output", {"type": "B_INT"})
    ET.SubElement(output_bint_elem, "Pressure").text = "true"
    ET.SubElement(output_bint_elem, "Velocity").text = "true"

    output_vint_elem = ET.SubElement(add_equation_elem, "Output", {"type": "V_INT"})
    ET.SubElement(output_vint_elem, "Pressure").text = "true"

    # LS block
    ls_elem = ET.SubElement(add_equation_elem, "LS", {"type": "NS"})
    lin_alg_elem = ET.SubElement(ls_elem, "Linear_algebra", {"type": "fsils"})
    ET.SubElement(lin_alg_elem, "Preconditioner").text = "fsils"

    ET.SubElement(ls_elem, "Max_iterations").text = "15"
    ET.SubElement(ls_elem, "NS_GM_max_iterations").text = "10"
    ET.SubElement(ls_elem, "NS_CG_max_iterations").text = "300"
    ET.SubElement(ls_elem, "Tolerance").text = "1e-3"
    ET.SubElement(ls_elem, "NS_GM_tolerance").text = "1e-3"
    ET.SubElement(ls_elem, "NS_CG_tolerance").text = "1e-3"
    ET.SubElement(ls_elem, "Absolute_tolerance").text = "1e-17"
    ET.SubElement(ls_elem, "Krylov_space_dimension").text = "250"

    # 7. Add_BC for each face
    for res, face in zip(resistances, cap_faces):
        bc_elem = ET.SubElement(add_equation_elem, "Add_BC", {"name": face})
        ET.SubElement(bc_elem, "Type").text = "Neu"
        ET.SubElement(bc_elem, "Time_dependence").text = "Resistance"
        ET.SubElement(bc_elem, "Value").text = str(res)
        
    for v, name in zip(inflows, inflow_names):
        bc_elem = ET.SubElement(add_equation_elem, "Add_BC", {"name": name})
        ET.SubElement(bc_elem, "Type").text = "Dir"
        ET.SubElement(bc_elem, "Time_dependence").text = "Steady"
        ET.SubElement(bc_elem, "Value").text = str(-1*v)
        ET.SubElement(bc_elem, "Profile").text = "Parabolic"
        ET.SubElement(bc_elem, "Impose_flux").text = "true"
        
    for fname in face_names:
        if fname == "walls_combined":
            bc_elem = ET.SubElement(add_equation_elem, "Add_BC", {"name": fname})
            # If it starts with cap_lpa or cap_rpa => Neu / Resistance
            if fname.startswith("wall"):
                ET.SubElement(bc_elem, "Type").text = "Dir"
                ET.SubElement(bc_elem, "Time_dependence").text = "Steady"
                ET.SubElement(bc_elem, "Value").text = "0.0"

    # 8. Pretty-print and write
    rough_string = ET.tostring(root, 'utf-8')
    reparsed = minidom.parseString(rough_string)
    pretty_xml = reparsed.toprettyxml(indent="  ")

    with open(output_xml_path, "w", encoding="utf-8") as f:
        f.write(pretty_xml)

    print(f"XML file written to: {output_xml_path}")

    # 9. Open result
    open_file_in_default_viewer(output_xml_path)


if __name__ == "__main__":
    project_directory = ""
    output_file = "solver.xml"
    areas, face_names = read_cap_info(project_directory)
    total_pressure_drop, left_ratio, inflows, inflow_names = prompt_BC(areas, face_names)
    resistances, cap_faces = calculate_resistance(total_pressure_drop, inflows, left_ratio, areas, face_names)
    create_mesh_xml(project_directory, output_file, resistances, cap_faces, inflows, inflow_names)
