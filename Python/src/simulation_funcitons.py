import os

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns 

from scipy.stats import t
import scipy.stats as stats
from scipy.spatial.distance import pdist, squareform
from scipy.optimize import curve_fit

import plotly.express as px
import plotly.graph_objects as go

from tqdm import tqdm

from src.helpers import print_difuso

class Simulador:
    def __init__(self, n, num_simulaciones, mu, variog_model, tau2, sigma2, phi, max_h, grid_type, min_x_grid, max_x_grid, min_y_grid, max_y_grid):
        """
        Ejemplo de uso:

        # Instanciamos la clase
        sim = Simulador(n=100,                     Número de puntos
                        num_simulaciones=300,       Numero de simulaciones a generar
                        mu=0,                       Vector de medias con este valor repetido
                        variog_model="gaussiano",   Tipo de modelo de variograma ('gaussiano', 'esferico', 'exponencial', 'sin_correlacion_espacial')
                        tau2=0,                     Efecto pepita
                        sigma2=5,                   Varianza
                        phi=30,                     Rango
                        max_h=100,                  Máxima distancia h para el gráfico de semivariograma
                        grid_type="random",         Tipo de enmallado ('random', 'regular')
                        min_x_grid=0,               Min x de la malla
                        max_x_grid=100,             Max x de la malla    
                        min_y_grid=0,               Min y de la malla
                        max_y_grid=100)             Max y de la malla
        sim.run()

        # Simulaciones (cada simulación es una fila):
        sim.samples
        """
        # Parámetros de la simulación
        self.n = n 
        self.num_simulaciones = num_simulaciones

        # Parámetros de la media
        self.mu = mu

        # Parámetros de la estructura de covarianza
        self.variog_model = variog_model
        self.tau2 = tau2
        self.sigma2 = sigma2
        self.phi = phi
        
        # Parámetros gráfico del modeo
        self.max_h = max_h

        # Parámetros del enmallado
        self.grid_type = grid_type
        self.min_x_grid = min_x_grid
        self.max_x_grid = max_x_grid
        self.min_y_grid = min_y_grid
        self.max_y_grid = max_y_grid
    

    def create_grid(self):
        """
        """
        if self.grid_type == "random":
            x_list = np.random.uniform(self.min_x_grid, self.max_x_grid, self.n)
            y_list = np.random.uniform(self.min_y_grid, self.max_y_grid, self.n)
            grid = pd.DataFrame({"x":x_list,
                                 "y":y_list})

        elif self.grid_type == "regular":

            # Verificamos que n sea un cuadrado perfecto
            if np.sqrt(self.n) != int(np.sqrt(self.n)):
                raise ValueError(f"Error: n ({self.n}) no es un cuadrado perfecto. Debe ser un cuadrado perfecto para usar un enmallado regular.")
    
            n_puntos = int(np.sqrt(self.n))

            x_list = np.linspace(self.min_x_grid, self.max_x_grid, n_puntos)
            y_list = np.linspace(self.min_y_grid, self.max_y_grid, n_puntos)
            
            x_list, y_list = np.meshgrid(x_list, y_list)
            
            grid = pd.DataFrame({"x":x_list.reshape(-1),
                                 "y":y_list.reshape(-1)})
            
        else:
            raise ValueError("Error: Introduzca un tipo permitico. El tipo debe ser 'random' o 'regular'")

        self.grid = grid

    def plot_grid(self):

        # Crear el gráfico de dispersión
        fig = px.scatter(x=self.grid.x, 
                         y=self.grid.y,
                         title="Gráfico de Grilla", 
                         labels={'x': 'Eje X', 'y': 'Eje Y'})

        # Actualizamos las dimensiones del gráfico
        fig.update_layout(
            width=800,  
            height=800, 
        )
        # Mostrar el gráfico
        fig.show()

    def semivariogram(self, h):
        """
        Calcula la semivarianza para un conjunto de distancias utilizando el modelo de semivariograma especificado.
        
        :param modelo: Tipo de modelo de semivariograma ('gaussiano', 'esferico', 'exponencial', 'sin_correlacion_espacial').
        :param distancias: Array de distancias entre los puntos.
        :param sigma2: Varianza estructural (σ^2).
        :param tau2: Varianza nugget (τ^2).
        :param phi: Rango del modelo (φ).
        :return: Array de semivarianzas para las distancias dadas.
        """
        if self.variog_model == 'gaussiano':
            return self.sigma2 * (1 - np.exp(-(h / self.phi) ** 2)) + self.tau2
        
        elif self.variog_model == 'esferico':
            # Aseguramos que la semivarianza no se haga decreciente
            semivarianza = np.zeros_like(h)
            
            # Para h <= phi, usamos la fórmula estándar
            semivarianza[h <= self.phi] = self.sigma2 * ((3 * h[h <= self.phi]) / (2 * self.phi) - (h[h <= self.phi] ** 3) / (2 * self.phi ** 3))
            
            # Para h > phi, la semivarianza es sigma^2 + tau^2
            semivarianza[h > self.phi] = self.sigma2 + self.tau2
            
            return semivarianza
        
        elif self.variog_model == 'exponencial':
            return self.sigma2 * (1 - np.exp(-h / self.phi)) + self.tau2
        
        elif self.variog_model == 'sin_correlacion_espacial':
            return np.full_like(h, self.tau2)  # Sin correlación espacial, todas las semivarianzas son tau2
        
        else:
            raise ValueError("Modelo de semivariograma no reconocido. Elige entre 'gaussiano', 'esferico', 'exponencial' o 'sin_correlacion_espacial'.")

    def generar_matriz_covarianzas(self):
        """
        Genera la matriz de covarianzas y varianzas para un conjunto de puntos usando un modelo de semivariograma.
        
        :param df: DataFrame con las coordenadas 'x' y 'y' de los puntos simulados.
        :param modelo: Tipo de modelo de semivariograma ('gaussiano', 'esferico', 'exponencial', 'sin_correlacion_espacial').
        :param sigma2: Varianza estructural (σ^2).
        :param tau2: Varianza nugget (τ^2).
        :param phi: Rango del modelo (φ).
        :return: Matriz de covarianzas de los puntos.
        """
        # Extraemos las coordenadas x, y del DataFrame
        coords = self.grid[['x', 'y']].values
        
        # Calculamos todas las distancias entre los puntos usando pdist y squareform
        distancias = pdist(coords)  # pdist calcula las distancias entre todos los pares de puntos
        
        # Calculamos la semivarianza para todas las distancias
        semivarianzas = self.semivariogram(h=distancias)
        
        if self.variog_model == "sin_correlacion_espacial":
            diagonal = np.full(self.n, self.sigma2)
            matriz_covarianzas = np.diag(diagonal)
        else:
            matriz_covarianzas = self.sigma2 - squareform(semivarianzas)  # Convertimos las semivarianzas a covarianzas
        
        # Devolvemos la matriz de covarianzas como un DataFrame
        self.cov_matrix = matriz_covarianzas

    def plot_cov_model(self):
        # Simulación de datos (ajusta según tus necesidades)
        h = np.linspace(0, self.max_h, 1000)
        gamma_h = self.semivariogram(h=h)

        # Crear el gráfico de dispersión
        plt.figure(figsize=(8, 8))  # Tamaño del gráfico

        # Graficar los datos
        plt.plot(h, gamma_h, color='blue', label="$\gamma(h)$")  # Puntos azules de dispersión

        # Crear el título con LaTeX
        plt.title(f"Modelo de Semivariograma {self.variog_model} con " + 
                f"$\\tau^2 = {self.tau2}$, " +
                f"$\\sigma^2 = {self.sigma2}$, " +
                f"$\\phi = {self.phi}$", fontsize=16)

        # Etiquetas de los ejes
        plt.xlabel(r"$h$", fontsize=14)
        plt.ylabel(r"$\gamma(h)$", fontsize=14)

        # Configuración del fondo y la cuadrícula
        plt.style.use('bmh')  # Fondo blanco
        plt.grid(True, which='both', axis='both', color='gray', linestyle='-', linewidth=0.5)

        # Añadimos el rango y la varianza
        plt.axhline(self.sigma2, color="orange", label="$\sigma^2$")
        plt.axvline(self.phi, color="green", label="$phi$")

        plt.legend()

        # Mostrar el gráfico
        plt.show()

    def simulate(self):
        # Definimos el vector de medias
        mu_vec = np.full(self.n, self.mu)
        # Generar las muestras
        self.samples = np.random.multivariate_normal(mu_vec, self.cov_matrix, self.num_simulaciones)

    def run(self):
        self.create_grid()
        #self.plot_grid()
        self.generar_matriz_covarianzas()
        #self.plot_cov_model()
        self.simulate()



class OneSampleTTestPower:
    def __init__(self, alpha, simulador, valores_grid_type, min_x_grid, max_x_grid, min_y_grid, max_y_grid, valores_mu, valores_variog_model, valores_sigma2, valores_tau2, valores_phi, valores_n, min_delta, max_delta, n_points_delta, n_bootstrap):
        
        # Parámetros de la simulación
        self.alpha = alpha
        self.min_delta = min_delta
        self.max_delta = max_delta
        self.n_points_delta = n_points_delta
        self.n_bootstrap = n_bootstrap

        # Simulador
        self.simulador = simulador

        # Variación de los parámetros
        self.valores_grid_type = valores_grid_type
        self.valores_mu = valores_mu
        self.valores_variog_model = valores_variog_model
        self.valores_sigma2 = valores_sigma2
        self.valores_tau2 = valores_tau2
        self.valores_phi = valores_phi
        self.valores_n = valores_n

        self.min_x_grid=min_x_grid,
        self.max_x_grid=max_x_grid,
        self.min_y_grid=min_y_grid,
        self.max_y_grid=max_y_grid,

    def get_power_function(self, n, mu, variog_model, tau2, sigma2, phi, grid_type):

        # Definimos el rango de los posibles delta's
        delta_range = np.linspace(self.min_delta, self.max_delta, self.n_points_delta)
        
        # Definimos una lista para guardar las potencias
        potencia = []

        for delta in tqdm(delta_range, desc=f"Simulaciones para n = {n}, mu = {mu}, variog_model = {variog_model}, tau2 = {tau2}, sigma2 = {sigma2}, phi = {phi}, grid_type = {grid_type}"):
            
            # Creamos un contador de rechazos
            rechazos = 0

            # Creamos un contador para la iteración bootstrap
            iter_bootstrap = 1

            for _ in range(self.n_bootstrap):
                sim = self.simulador(n=n,
                                    num_simulaciones=1,
                                    mu=mu + delta,
                                    variog_model=variog_model,
                                    tau2=tau2, 
                                    sigma2=sigma2, 
                                    phi=phi,
                                    max_h=100, 
                                    grid_type=grid_type, 
                                    min_x_grid=self.min_x_grid, 
                                    max_x_grid=self.max_x_grid, 
                                    min_y_grid=self.min_y_grid, 
                                    max_y_grid=self.max_y_grid)
                sim.run()
                # Hacemos la prueba
                t_stat, p_value = stats.ttest_1samp(sim.samples.reshape(-1), mu, alternative="greater")

                # Vemos si se rechaza la hipótesis nula
                if p_value < self.alpha:
                    rechazos += 1

                # Imprimimos la iteración que vamos en el bootstrap
                # print_difuso(f"{iter_bootstrap} / {self.n_bootstrap}")
                # iter_bootstrap = iter_bootstrap + 1
                # Comentamos estas lineas porque consumen mucho tiempo

            # Añadimos la potencia al vector de potencias
            potencia.append(rechazos / self.n_bootstrap)

        return delta_range, potencia

    def run_power_functions(self):
        resultados = []

        for n in self.valores_n:
            for mu in self.valores_mu:
                for variog_model in self.valores_variog_model:
                    for tau2 in self.valores_tau2:
                        for sigma2 in self.valores_sigma2:
                            for phi in self.valores_phi:
                                for grid_type in self.valores_grid_type:
                                    delta_range, potencia = self.get_power_function(n, mu, variog_model, tau2, sigma2, phi, grid_type)

                                    dic_resultados = {"alpha": self.alpha,
                                                      "n": n,
                                                      "mu": mu,
                                                      "variog_model": variog_model,
                                                      "tau2": tau2,
                                                      "sigma2": sigma2,
                                                      "phi": phi,
                                                      "grid_type" : grid_type,
                                                      "min_x_grid" : self.min_x_grid,
                                                      "max_x_grid" : self.max_x_grid,
                                                      "min_y_grid" : self.min_y_grid,
                                                      "max_y_grid" : self.max_y_grid,
                                                      "min_delta" : self.min_delta,
                                                      "max_delta" : self.max_delta,
                                                      "n_points_delta" : self.n_points_delta,
                                                      "n_bootstrap" : self.n_bootstrap,
                                                      "delta_range" : delta_range,
                                                      "potencia" : potencia}
                                    resultados.append(dic_resultados)

        self.resultados_totales = pd.DataFrame(resultados)
        self.resultados_totales.to_csv(f"resultados_totales_bootstrap_{self.n_bootstrap}.csv", index=False)

    def run(self):
        self.run_power_functions()