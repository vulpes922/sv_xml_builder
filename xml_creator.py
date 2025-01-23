import os
import glob
import xml.etree.ElementTree as ET
from xml.dom import minidom
import sys
import subprocess

def open_file_in_default_viewer(file_path):
    """
    Open the specified file in the system's default application.
    Works on Windows, macOS, and Linux (xdg-open).
    """
    if sys.platform.startswith("win"):
        os.startfile(file_path)
    elif sys.platform == "darwin":
        subprocess.run(["open", file_path])
    else:
        subprocess.run(["xdg-open", file_path])

def create_mesh_xml(project_dir, output_xml_path):
    """
    Create an XML file with:
      <svMultiPhysicsFile version="0.1">
        <GeneralSimulationParameters>...</GeneralSimulationParameters>
        <Add_mesh name="msh">...</Add_mesh>
      </svMultiPhysicsFile>
    
    and then open it in the default viewer.
    """

    # 1. Locate the mesh and mesh-surfaces folders
    mesh_dir = os.path.join(project_dir, "mesh")
    mesh_surfaces_dir = os.path.join(mesh_dir, "mesh-surfaces")

    if not os.path.isdir(mesh_dir):
        raise FileNotFoundError(f"Mesh folder not found at: {mesh_dir}")
    if not os.path.isdir(mesh_surfaces_dir):
        raise FileNotFoundError(f"Mesh-surfaces folder not found at: {mesh_surfaces_dir}")

    # 2. Gather all .vtp files
    vtp_files = glob.glob(os.path.join(mesh_surfaces_dir, "*.vtp"))
    if not vtp_files:
        print(f"No .vtp files found in: {mesh_surfaces_dir}")

    # 3. Build the root element: <svMultiPhysicsFile>
    root = ET.Element("svMultiPhysicsFile", {"version": "0.1"})
    root.tail = " "

    # -------------------------------------------------------------------------
    # 4. Add the <GeneralSimulationParameters> section
    # -------------------------------------------------------------------------
    gsp = ET.SubElement(root, "GeneralSimulationParameters")
    gsp.tail = " "

    # Populate each sub-tag
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

    # -------------------------------------------------------------------------
    # 5. Add the <Add_mesh> element
    # -------------------------------------------------------------------------
    add_mesh_elem = ET.SubElement(root, "Add_mesh", {"name": "msh"})
    add_mesh_elem.tail = " "

    # The primary volume mesh:
    mesh_file_path = ET.SubElement(add_mesh_elem, "Mesh_file_path")
    mesh_file_path.text = "mesh/mesh-complete.mesh.vtu"
    mesh_file_path.tail = " "
    

    # For each .vtp file, create <Add_face>
    for vtp_file in vtp_files:
        # face_name = filename (without extension)
        face_name = os.path.splitext(os.path.basename(vtp_file))[0]

        add_face_elem = ET.SubElement(add_mesh_elem, "Add_face", {"name": face_name})
        face_file_path_elem = ET.SubElement(add_face_elem, "Face_file_path")
        face_file_path_elem.text = f"mesh/mesh-surfaces/{os.path.basename(vtp_file)}"
        add_face_elem.tail = " "
    
    add_equation_elem = ET.SubElement(root, "Add_equation", {"type": "fluid"})

    # 1) Simple sub-elements (Coupled, Min_iterations, etc.)
    ET.SubElement(add_equation_elem, "Coupled").text = "true"
    ET.SubElement(add_equation_elem, "Min_iterations").text = "3"
    ET.SubElement(add_equation_elem, "Max_iterations").text = "5"
    ET.SubElement(add_equation_elem, "Tolerance").text = "1e-11"
    ET.SubElement(add_equation_elem, "Backflow_stabilization_coefficient").text = "0.2"

    # 2) Fluid properties
    ET.SubElement(add_equation_elem, "Density").text = "1.06"

    # Viscosity element with 'model="Constant"' attribute
    viscosity_elem = ET.SubElement(add_equation_elem, "Viscosity", {"model": "Constant"})
    ET.SubElement(viscosity_elem, "Value").text = "0.04"

    # 3) Output sections
    #    <Output type="Spatial"> ... </Output>
    output_spatial_elem = ET.SubElement(add_equation_elem, "Output", {"type": "Spatial"})
    ET.SubElement(output_spatial_elem, "Velocity").text = "true"
    ET.SubElement(output_spatial_elem, "Pressure").text = "true"
    ET.SubElement(output_spatial_elem, "Traction").text = "true"
    ET.SubElement(output_spatial_elem, "Vorticity").text = "true"
    ET.SubElement(output_spatial_elem, "Divergence").text = "true"
    ET.SubElement(output_spatial_elem, "WSS").text = "true"

    #    <Output type="B_INT"> ... </Output>
    output_bint_elem = ET.SubElement(add_equation_elem, "Output", {"type": "B_INT"})
    ET.SubElement(output_bint_elem, "Pressure").text = "true"
    ET.SubElement(output_bint_elem, "Velocity").text = "true"
    #    <Output type="V_INT"> ... </Output>
    output_vint_elem = ET.SubElement(add_equation_elem, "Output", {"type": "V_INT"})
    ET.SubElement(output_vint_elem, "Pressure").text = "true"
    # 4) <LS type="NS"> block
    ls_elem = ET.SubElement(add_equation_elem, "LS", {"type": "NS"})

    #    <Linear_algebra type="fsils">
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

    face_names = []

    for vtp_file in vtp_files:
        face_name = os.path.splitext(os.path.basename(vtp_file))[0]
        face_names.append(face_name)
    
    for fname in face_names:
        # e.g. if you want this BC *only* for faces that have the word "wall"
        # if "wall" not in fname.lower():
        #     continue

        bc_elem = ET.SubElement(add_equation_elem, "Add_BC", {"name": fname})
        ET.SubElement(bc_elem, "Type").text = "Dir"
        ET.SubElement(bc_elem, "Time_dependence").text = "Steady"
        ET.SubElement(bc_elem, "Value").text = "0.0"
        bc_elem.tail = " "
    
    # -------------------------------------------------------------------------
    # 6. Pretty-print and write to file
    # -------------------------------------------------------------------------
    rough_string = ET.tostring(root, 'utf-8')
    reparsed = minidom.parseString(rough_string)
    pretty_xml = reparsed.toprettyxml(indent="  ")

    with open(output_xml_path, "w", encoding="utf-8") as f:
        f.write(pretty_xml)

    print(f"XML file written to: {output_xml_path}")

    # -------------------------------------------------------------------------
    # 7. Open the resulting XML file in the default viewer
    # -------------------------------------------------------------------------
    open_file_in_default_viewer(output_xml_path)


if __name__ == "__main__":
    # Example usage:
    # Replace this path with the actual project directory on your system
    project_directory = "//wsl.localhost/Ubuntu/home/bryan922/svMultiPhysics/Projects/fontan/250122_example_fontan"
    
    # Output XML file path
    output_file = "solver.xml"
    
    # Create the XML
    create_mesh_xml(project_directory, output_file)
