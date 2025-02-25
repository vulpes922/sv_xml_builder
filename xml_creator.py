import os
import glob
import xml.etree.ElementTree as ET
from xml.dom import minidom
import sys
import subprocess
import tkinter as tk
from tkinter import filedialog, messagebox
import platform

def open_file_in_default_viewer(file_path):
    # 절대 경로로 변환
    file_path = os.path.abspath(file_path)
    # 경로 형식 통일 (선택사항)
    if platform.system() == 'Windows':
        path_sep = '\\'
    else:
        path_sep = '/'
    file_path = file_path.replace("/", "\\")   # there is a possibility that this depends on os linux, windows or mac
    file_path = file_path.replace("\\", path_sep) 
    if not os.path.exists(file_path):
        print(f"파일을 찾을 수 없습니다: {file_path}")
        return
    if sys.platform.startswith("win"):
        os.startfile(file_path)
    elif sys.platform == "darwin":
        subprocess.run(["open", file_path])
    else:
        subprocess.run(["xdg-open", file_path])



def read_cap_info(project_dir):
    """
    Read 'cap_info.txt' from the project directory.
    Returns two parallel lists: areas (floats) and face_names (strings).
    """
    cap_file_path = os.path.join(project_dir, "cap_info.txt")
    if not os.path.isfile(cap_file_path):
        raise FileNotFoundError(f"Could not find 'cap_info.txt' in: {project_dir}")
    data = []
    with open(cap_file_path, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('#'):
                continue
            parts = line.split()
            if len(parts) != 2:
                raise ValueError(f"Malformed line in cap_info.txt: '{line}'")
            face_name, area_str = parts
            area = float(area_str)
            data.append((face_name, area))
    # sort alphabetically by face_name
    data.sort(key=lambda x: x[0])
    face_names = [item[0] for item in data]
    areas = [item[1] for item in data]
    return areas, face_names

def calculate_resistance(total_pressure_drop, inflows, left_ratio, areas, face_names):
    """
    Calculate resistances for cap faces given the total pressure drop,
    inflows (for non-lpa/rpa faces), and the left/right pulmonary flow ratio.
    """
    total_inflow = sum(inflows)
    lpa = total_inflow * left_ratio/100
    rpa = total_inflow - lpa

    pressure = total_pressure_drop * 1333.2  # convert mmHg to dynes/cm^2

    # Compute total areas for left and right cap faces
    left_total_area = sum(a for a, face in zip(areas, face_names) if face.startswith("cap_lpa"))
    right_total_area = sum(a for a, face in zip(areas, face_names) if face.startswith("cap_rpa"))

    area_ratios = []
    cap_faces = []
    for a, face in zip(areas, face_names):
        if face.startswith("cap_lpa"):
            area_ratios.append(a / left_total_area)
            cap_faces.append(face)
        elif face.startswith("cap_rpa"):
            area_ratios.append(a / right_total_area)
            cap_faces.append(face)

    resistances = []
    for a, face in zip(area_ratios, cap_faces):
        if face.startswith("cap_lpa"):
            resistances.append(pressure / (a * lpa))
        elif face.startswith("cap_rpa"):
            resistances.append(pressure / (a * rpa))
    return resistances, cap_faces

def create_mesh_xml(project_dir, output_xml_path, resistances, cap_faces, inflows, inflow_names):
    """
    Create the XML file for the solver.
    It uses the mesh files in the 'mesh' and 'mesh-surfaces' subfolders of project_dir.
    """
    mesh_dir = os.path.join(project_dir, "mesh")
    mesh_surfaces_dir = os.path.join(mesh_dir, "mesh-surfaces")
    if not os.path.isdir(mesh_dir):
        raise FileNotFoundError(f"Mesh folder not found at: {mesh_dir}")
    if not os.path.isdir(mesh_surfaces_dir):
        raise FileNotFoundError(f"Mesh-surfaces folder not found at: {mesh_surfaces_dir}")

    # Get .vtp files that start with 'cap' and also include walls_combined.vtp if present.
    all_vtp_files = glob.glob(os.path.join(mesh_surfaces_dir, "*.vtp"))
    filtered_vtp_files = [vtp for vtp in all_vtp_files if os.path.basename(vtp).startswith("cap")]
    walls_combined_path = os.path.join(mesh_dir, "walls_combined.vtp")
    if os.path.isfile(walls_combined_path):
        filtered_vtp_files.append(walls_combined_path)
    filtered_vtp_files.sort(key=lambda p: os.path.basename(p).lower())

    root = ET.Element("svMultiPhysicsFile", {"version": "0.1"})
    root.tail = " "

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

    add_mesh_elem = ET.SubElement(root, "Add_mesh", {"name": "msh"})
    add_mesh_elem.tail = " "
    mesh_file_path = ET.SubElement(add_mesh_elem, "Mesh_file_path")
    mesh_file_path.text = "mesh/mesh-complete.mesh.vtu"
    mesh_file_path.tail = " "

    face_names = []
    for vtp_file in filtered_vtp_files:
        base_name = os.path.splitext(os.path.basename(vtp_file))[0]
        add_face_elem = ET.SubElement(add_mesh_elem, "Add_face", {"name": base_name})
        face_file_path_elem = ET.SubElement(add_face_elem, "Face_file_path")
        if base_name != "walls_combined":
            face_file_path_elem.text = f"mesh/mesh-surfaces/{os.path.basename(vtp_file)}"
        else:
            face_file_path_elem.text = f"mesh/{os.path.basename(vtp_file)}"
        add_face_elem.tail = " "
        face_names.append(base_name)

    add_equation_elem = ET.SubElement(root, "Add_equation", {"type": "fluid"})
    ET.SubElement(add_equation_elem, "Coupled").text = "true"
    ET.SubElement(add_equation_elem, "Min_iterations").text = "3"
    ET.SubElement(add_equation_elem, "Max_iterations").text = "5"
    ET.SubElement(add_equation_elem, "Tolerance").text = "1e-11"
    ET.SubElement(add_equation_elem, "Backflow_stabilization_coefficient").text = "0.2"
    ET.SubElement(add_equation_elem, "Density").text = "1.06"
    viscosity_elem = ET.SubElement(add_equation_elem, "Viscosity", {"model": "Constant"})
    ET.SubElement(viscosity_elem, "Value").text = "0.04"

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

    for res, face in zip(resistances, cap_faces):
        bc_elem = ET.SubElement(add_equation_elem, "Add_BC", {"name": face})
        ET.SubElement(bc_elem, "Type").text = "Neu"
        ET.SubElement(bc_elem, "Time_dependence").text = "Resistance"
        ET.SubElement(bc_elem, "Value").text = str(res)

    for v, name in zip(inflows, inflow_names):
        bc_elem = ET.SubElement(add_equation_elem, "Add_BC", {"name": name})
        ET.SubElement(bc_elem, "Type").text = "Dir"
        ET.SubElement(bc_elem, "Time_dependence").text = "Steady"
        ET.SubElement(bc_elem, "Value").text = str(-1 * v)
        ET.SubElement(bc_elem, "Profile").text = "Parabolic"
        ET.SubElement(bc_elem, "Impose_flux").text = "true"

    for fname in face_names:
        if fname == "walls_combined":
            bc_elem = ET.SubElement(add_equation_elem, "Add_BC", {"name": fname})
            if fname.startswith("wall"):
                ET.SubElement(bc_elem, "Type").text = "Dir"
                ET.SubElement(bc_elem, "Time_dependence").text = "Steady"
                ET.SubElement(bc_elem, "Value").text = "0.0"

    rough_string = ET.tostring(root, 'utf-8')
    reparsed = minidom.parseString(rough_string)
    pretty_xml = reparsed.toprettyxml(indent="  ")

    with open(output_xml_path, "w", encoding="utf-8") as f:
        f.write(pretty_xml)

    print(f"XML file written to: {output_xml_path}")
    open_file_in_default_viewer(output_xml_path)

class SolverGUI(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("Solver XML Generator")
        self.geometry("600x400")
        self.project_directory = None
        self.areas = []
        self.face_names = []
        self.inflow_faces = []
        self.inflow_entries = {}

        # Create UI elements
        self.create_widgets()

    def create_widgets(self):
        # --- Project Directory Selection ---
        dir_frame = tk.Frame(self)
        dir_frame.pack(pady=10)
        tk.Label(dir_frame, text="Project Directory:").pack(side=tk.LEFT)
        self.dir_label = tk.Label(dir_frame, text="None selected", fg="blue")
        self.dir_label.pack(side=tk.LEFT, padx=5)
        tk.Button(dir_frame, text="Browse", command=self.browse_directory).pack(side=tk.LEFT)

        # --- Total Pressure Drop and Left Ratio Inputs ---
        input_frame = tk.Frame(self)
        input_frame.pack(pady=10)
        tk.Label(input_frame, text="Total Pressure Drop (mmHg):").grid(row=0, column=0, sticky="e")
        self.total_pressure_drop_entry = tk.Entry(input_frame)
        self.total_pressure_drop_entry.grid(row=0, column=1, padx=5, pady=5)
        tk.Label(input_frame, text="Pulmonary Flow Ratio (percentage to left):").grid(row=1, column=0, sticky="e")
        self.left_ratio_entry = tk.Entry(input_frame)
        self.left_ratio_entry.grid(row=1, column=1, padx=5, pady=5)

        # --- Inflow Velocities ---
        self.inflow_frame = tk.LabelFrame(self, text="Inflow Velocities (for non cap_lpa/cap_rpa faces)")
        self.inflow_frame.pack(fill="both", expand=True, padx=10, pady=10)

        # --- Generate XML Button ---
        tk.Button(self, text="Generate XML", command=self.generate_xml).pack(pady=10)

    def browse_directory(self):
        directory = filedialog.askdirectory()
        if directory:
            self.project_directory = directory
            self.dir_label.config(text=directory)
            try:
                self.areas, self.face_names = read_cap_info(self.project_directory)
                # Clear any existing entries
                for widget in self.inflow_frame.winfo_children():
                    widget.destroy()
                self.inflow_faces = []
                self.inflow_entries = {}
                row = 0
                # For every face that is not cap_lpa or cap_rpa, create an entry
                for face, area in zip(self.face_names, self.areas):
                    if not (face.startswith("cap_lpa") or face.startswith("cap_rpa")):
                        self.inflow_faces.append(face)
                        tk.Label(self.inflow_frame, text=f"{face} (area = {area}):").grid(row=row, column=0, sticky="e", padx=5, pady=2)
                        entry = tk.Entry(self.inflow_frame)
                        entry.grid(row=row, column=1, padx=5, pady=2)
                        self.inflow_entries[face] = entry
                        row += 1
            except Exception as e:
                messagebox.showerror("Error", str(e))

    def generate_xml(self):
        try:
            total_pressure_drop = float(self.total_pressure_drop_entry.get())
            left_ratio = float(self.left_ratio_entry.get())
        except ValueError:
            messagebox.showerror("Input error", "Please enter valid numeric values for pressure drop and ratio.")
            return

        inflows = []
        inflow_names = []
        for face in self.inflow_faces:
            try:
                velocity = float(self.inflow_entries[face].get())
            except ValueError:
                messagebox.showerror("Input error", f"Please enter a valid numeric value for inflow velocity for face {face}.")
                return
            inflows.append(velocity)
            inflow_names.append(face)
        try:
            resistances, cap_faces = calculate_resistance(total_pressure_drop, inflows, left_ratio, self.areas, self.face_names)
            output_file = os.path.join(self.project_directory, "solver.xml")
            create_mesh_xml(self.project_directory, output_file, resistances, cap_faces, inflows, inflow_names)
            messagebox.showinfo("Success", f"XML file generated at: {output_file}")
        except Exception as e:
            messagebox.showerror("Error", str(e))

if __name__ == "__main__":
    app = SolverGUI()
    app.mainloop()
