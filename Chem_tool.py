# -*- coding: utf-8 -*-
"""
Created on Fri Oct 22 09:05:18 2026

@author: Carlos F. Marcos
"""

import streamlit as st

# --- 1. CONFIGURACIÓN DE LA PÁGINA ---
st.set_page_config(
    page_title="Cuaderno y Esquemas de Laboratorio",
    layout="wide",
    page_icon="🧪",
    initial_sidebar_state="expanded" 
)

# --- 2. IMPORTACIÓN DE MÓDULOS ---
try:
    from Reaction_Scheme_Exporter import render_reaction_scheme_exporter as render_reaction_scheme_app
    from Reaction_Lab_Notebook import render_reaction_lab_notebook as render_labnotebook_app
except ImportError as e:
    st.error(f"Error cargando módulos de la suite: {e}")

# --- 3. BARRA LATERAL (SIDEBAR) ---
with st.sidebar:
    st.title("🧪 Suite de Reacciones")
    st.markdown("Herramientas de Química y Estequiometría")
    st.divider()
    
    # Selector de herramientas en español
    opcion = st.radio(
        "Selecciona una herramienta:", 
        [
            "Figuras y Esquemas de Reacción", 
            "Libreta de Laboratorio"
        ]
    )

# --- Colores dinámicos para el sidebar según la herramienta seleccionada ---
tool_colors = {
    "Figuras y Esquemas de Reacción": "#e3f2fd",  # Azul claro suave
    "Libreta de Laboratorio": "#f1f8e9"          # Verde claro suave
}

def set_sidebar_color(option):
    color = tool_colors.get(option, "#ffffff")
    st.markdown(f"""
        <style>
            [data-testid="stSidebar"] {{
                background-color: {color};
            }}
        </style>
    """, unsafe_allow_html=True)

set_sidebar_color(opcion)

# --- 4. ENRUTAMIENTO (ROUTING) ---
if opcion == "Figuras y Esquemas de Reacción":
    render_reaction_scheme_app()
elif opcion == "Libreta de Laboratorio":
    render_labnotebook_app()