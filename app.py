import numpy as np
import plotly.graph_objects as go
import streamlit as st

# --- CONFIGURACIÓN DE LA PÁGINA ---
st.set_page_config(
    page_title="Diagrama de Ostwald", page_icon="🔥", layout="wide"
)

st.title("🔥 Generador de Diagrama de Ostwald")

# --- PANEL LATERAL (CONTROLES Y PARÁMETROS) ---
st.sidebar.header("1. Parámetros del Combustible")
combustible_nombre = st.sidebar.text_input(
    "Nombre del Combustible", value="Metano (CH4)"
)
n = st.sidebar.number_input(
    "Átomos de Carbono (n)", min_value=0.1, value=1.0, step=0.1
)
m = st.sidebar.number_input(
    "Átomos de Hidrógeno (m)", min_value=0.0, value=4.0, step=0.1
)

st.sidebar.markdown("---")
st.sidebar.header("2. Lectura de Chimenea")
activar_medicion = st.sidebar.checkbox("Ingresar punto medido", value=True)

if activar_medicion:
  o2_medido = st.sidebar.number_input(
      "O₂ medido [omega]",
      min_value=0.0,
      max_value=21.0,
      value=3.0,
      step=0.1,
  )
  co2_medido = st.sidebar.number_input(
      "CO₂ medido [alpha]",
      min_value=0.0,
      max_value=20.0,
      value=8.0,
      step=0.1,
  )
else:
  o2_medido, co2_medido = None, None

# --- PARÁMETROS DE OSTWALD EN EL PANEL LATERAL ---
st.sidebar.markdown("---")
st.sidebar.header("3. Parámetros de Ostwald")

# 1. Parámetros fundamentales
alpha_m = n / (n + 3.762 * (n + m / 4))

if activar_medicion and o2_medido is not None and co2_medido is not None:
  omega_val = o2_medido / 100.0
  alpha_val = co2_medido / 100.0

  # Cálculo exacto de n_s según las fórmulas de cátedra
  denominador_ns = 1 - 4.762 * omega_val - (3.762 / 2) * alpha_val

  if denominador_ns > 0:
    num_ns = n + 3.762 * (n / 2 + m / 4)
    ns_val = num_ns / denominador_ns
    x_val = n - alpha_val * ns_val
    z_val = omega_val * ns_val
    gamma_val = x_val / ns_val
  else:
    ns_val, x_val, z_val, gamma_val = 0.0, 0.0, 0.0, 0.0
else:
  omega_val, alpha_val = 0.0, 0.0
  ns_val, x_val, z_val, gamma_val = 0.0, 0.0, 0.0, 0.0

st.sidebar.markdown(f"**$\\alpha$ =** `{alpha_val:.4f}`")
st.sidebar.markdown(f"**$\\alpha_m$ =** `{alpha_m:.4f}`")
st.sidebar.markdown(f"**$\\omega$ =** `{omega_val:.4f}`")
st.sidebar.markdown(f"**$n_s$ (moles humos secos)=** `{ns_val:.4f}`")
st.sidebar.markdown(f"**$x$ (moles CO)=** `{x_val:.4f}`")
st.sidebar.markdown(f"**$z$ (moles O₂)=** `{z_val:.4f}`")
st.sidebar.markdown(f"**$\\gamma$ =** `{gamma_val:.4f}`")


# --- SECCIÓN DE ECUACIÓN GENERAL DE OSTWALD ---
st.subheader("📝 Ecuación General de Reacción de Ostwald")
st.latex(
    r"C_n H_m + \left(n + \frac{m}{4} - \frac{x}{2} + z\right)O_2 + 3.762\left(n"
    r" + \frac{m}{4} - \frac{x}{2} + z\right)N_2 ="
    r"(n-x)CO_2 + x CO + \frac{m}{2}H_2O + z O_2 + 3.762\left(n + \frac{m}{4} -"
    r" \frac{x}{2} + z\right)N_2"
)

st.markdown("---")

# --- CONSTRUCCIÓN DEL GRÁFICO ---
omega = np.linspace(0, 0.21, 500)
fig = go.Figure()

# Recta de Grebel (CO = 0%)
alpha_grebel = alpha_m * (1 - 4.762 * omega)
fig.add_trace(
    go.Scatter(
        x=omega,
        y=alpha_grebel,
        mode="lines",
        name="Recta de Grebel (CO = 0%)",
        line=dict(color="black", width=3),
        hovertemplate=(
            "<b>Recta de Grebel</b><br>O₂: %{x:.4f}<br>CO₂: %{y:.4f}<extra></extra>"
        ),
    )
)

# Líneas de CO constante
valores_gamma = [i / 100.0 for i in range(1, 15)]
for g in valores_gamma:
  alpha_co = (3.762 / 2 * alpha_m - 1) * g + alpha_m * (1 - 4.762 * omega)
  mask = (alpha_co >= 0) & (alpha_co <= alpha_grebel)
  if np.any(mask):
    fig.add_trace(
        go.Scatter(
            x=omega[mask],
            y=alpha_co[mask],
            mode="lines",
            name=f"CO = {int(g*100)}%",
            line=dict(color="rgba(0, 102, 204, 0.35)", width=1.2, dash="dash"),
            hovertemplate=(
                f"<b>CO = {int(g*100)}%</b><br>O₂: %{{x:.4f}}<br>CO₂:"
                " %{y:.4f}<extra></extra>"
            ),
        )
    )

# Líneas de exceso de aire (e)
polo_x, polo_y = 1.0, -2.0
valores_e = [
    -0.2,
    -0.15,
    -0.1,
    -0.05,
    0,
    0.05,
    0.1,
    0.15,
    0.2,
    0.3,
    0.5,
    1.0,
    2.0,
    5.0,
    10.0,
]

for e in valores_e:
  num_w = e * (n + m / 4)
  den_w = n + 3.762 * (n + m / 4) * (1 + e)
  omega_g = num_w / den_w
  alpha_g = alpha_m * (1 - 4.762 * omega_g)

  m_e = (alpha_g - polo_y) / (omega_g - polo_x)
  b_e = polo_y - m_e * polo_x

  a_fit = m_e * omega + b_e
  mask_e = (a_fit >= 0) & (a_fit <= alpha_grebel + 0.0005)

  if np.any(mask_e):
    lbl_e = f"e = {e}" if e != 0 else "e = 0 (Estequiométrico)"
    fig.add_trace(
        go.Scatter(
            x=omega[mask_e],
            y=a_fit[mask_e],
            mode="lines",
            name=lbl_e,
            line=dict(color="rgba(34, 139, 34, 0.75)", width=1.3),
            hovertemplate=(
                f"<b>{lbl_e}</b><br>O₂: %{{x:.4f}}<br>CO₂:"
                " %{y:.4f}<extra></extra>"
            ),
        )
    )

# Graficar Punto Medido
if activar_medicion and o2_medido is not None and co2_medido is not None:
  fig.add_trace(
      go.Scatter(
          x=[omega_val],
          y=[alpha_val],
          mode="markers",
          name="Punto Medido",
          marker=dict(color="red", size=12, symbol="cross", line=dict(width=2)),
          hovertemplate=(
              f"<b>Punto Medido</b><br>O₂: {o2_medido}%<br>CO₂:"
              f" {co2_medido}%<extra></extra>"
          ),
      )
  )

# Estética del gráfico
fig.update_layout(
    title=(
        "<b>Diagrama de Ostwald Interactivo -"
        f" {combustible_nombre}</b>"
    ),
    xaxis_title="Fracción molar de O₂ en humos secos (ω)",
    yaxis_title="Fracción molar de CO₂ en humos secos (α)",
    xaxis=dict(range=[0, 0.215], showgrid=True, gridcolor="#E5E5E5"),
    yaxis=dict(range=[0, alpha_m * 1.08], showgrid=True, gridcolor="#E5E5E5"),
    hovermode="closest",
    template="plotly_white",
    height=650,
    showlegend=False,
)

# --- MOSTRAR RESULTADOS EN LA INTERFAZ ---
col1, col2 = st.columns([3, 1])

with col1:
  st.plotly_chart(fig, use_container_width=True)

with col2:
  st.subheader("📊 Diagnóstico de Chimenea")
  

  if activar_medicion:
    st.markdown("---")
    

    a_grebel_punto = alpha_m * (1 - 4.762 * omega_val)
    if alpha_val > a_grebel_punto + 0.0005:
      st.error(
          "⚠️ **Punto Imposible:** Las concentraciones están por encima de la"
          " recta de combustión completa."
      )
    else:
      st.metric(
          label="CO en humos (γ)",
          value=f"{gamma_val*100:.2f} %",
          delta=f"{x_val:.3f} moles/mol comb",
      )
      e_calc = (
          ((z_val - (x_val/2))/(n+(m/4))) if (0.21 - omega_val) > 0 else 0.0
      )
      st.metric(
          label="Exceso de Aire Aprox. (e)", value=f"{e_calc*100:.1f} %"
      )

      if gamma_val > 0.0001:
        st.warning("⚠️ **Combustión Incompleta:** Se detectó CO.")
      else:
        st.success("✅ **Combustión Completa:** Prácticamente 0% de CO.")