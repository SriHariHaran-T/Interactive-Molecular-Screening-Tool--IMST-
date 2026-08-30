# Interactive Molecular Screening Tool (IMST) - Pitch & User Manual

Welcome to the **Interactive Molecular Screening Tool (IMST)**! 

IMST is a computational drug-screening web application designed to streamline the early stages of Computer-Aided Drug Design (CADD). By focusing on ADMET-based oral drug screening, IMST instantly predicts whether a molecule is likely to be effectively absorbed by the human gut and whether it fits a general oral drug profile.

---

## 🚀 How to Access IMST

You can access the live interactive web application right now from your browser. No installation or setup is required!

**👉 Click here to launch the app: [https://imstproject.streamlit.app/](https://imstproject.streamlit.app/)**

---

## 🌟 Advantages: Why Choose IMST?

- **Instant Insights**: Get immediate predictions on oral bioavailability and drug-likeness without running complex local software.
- **Scientifically Grounded**: Built on established cheminformatics principles, specifically Lipinski’s Rule of Five and Veber’s Rule.
- **Highly Interactive**: Move beyond static data. Visualize your molecules in full 3D, right in your browser.
- **User-Friendly**: Designed with an intuitive interface that accommodates both beginners learning about drug design and researchers needing quick validations.
- **Accessible Anywhere**: As a cloud-hosted web app, you can screen molecules from any device, anytime.

---

## 🎯 Key Uses

IMST is designed to be a versatile tool for various users in the scientific community:

1. **Early-Stage Drug Discovery**: Quickly filter out unviable molecular candidates before investing in expensive and time-consuming in-vitro or in-vivo testing.
2. **Educational Tool**: An excellent resource for chemistry, biology, and data science students to visualize and understand molecular properties and pharmacokinetics.
3. **Rapid Prototyping**: Researchers can rapidly test custom SMILES strings to see how structural modifications impact a molecule's predicted absorption.

---

## 🛠️ Core Features

- **ADMET Screening Engine**: Automatically computes core physicochemical descriptors (Molecular Weight, LogP, TPSA, H-bond donors/acceptors, rotatable bonds).
- **Rule-Based Evaluation**: 
  - *Lipinski’s Rule of Five*: Evaluates general oral drug-likeness.
  - *Veber’s Rule*: Predicts human gut absorption based on TPSA and rotatable bonds.
- **Live Verdict Badges**: Instantly see easy-to-read badges indicating "Likely absorbed by the gut" or "Fits general oral drug profile".
- **Interactive 3D Viewer**: A rotatable, zoomable 3D molecular structure viewer powered by 3Dmol.js.
- **Molecular Docking Pipeline**: Perform molecular docking against protein targets (e.g., HIV-1 Protease) using AutoDock Vina, complete with protein-ligand contact analysis and 3D pose visualization.

---

## 📖 How to Use the App

Using IMST is straightforward and requires just a few clicks:

### Method 1: Explore Known Drugs
Perfect for getting a feel for the tool or studying well-known compounds.
1. Open the [IMST Web App](https://imstproject.streamlit.app/).
2. Under "Choose input method", select **Pick from list**.
3. Use the dropdown menu to select from our curated starter dataset of 15 verified real drugs (e.g., Aspirin, Ibuprofen, Caffeine).
4. Watch as the property table, verdict badges, and 3D structure automatically update!

### Method 2: Test Custom Molecules
Ideal for testing new ideas and novel compounds.
1. Open the [IMST Web App](https://imstproject.streamlit.app/).
2. Under "Choose input method", select **Enter my own SMILES**.
3. Paste or type the SMILES string of your molecule into the input field.
4. Hit enter, and IMST will instantly compute the descriptors, render the 3D structure, and deliver the ADMET verdicts.

### Method 3: Molecular Docking
Screen your molecule against a protein target.
1. Navigate to the **Docking** section in the web app.
2. Ensure your SMILES string is loaded.
3. Initiate the docking process (powered by AutoDock Vina).
4. Analyze the resulting docking score, binding pose, and protein-ligand interactions in the 3D viewer.

---

*IMST is actively being developed by the First-Year Integrated M.Sc. Data Science program at TCE Madurai. It now includes full molecular docking capabilities!*
