import yaml
import time
import pyfiglet


def load_config(ruta_archivo):
    """Lee un archivo YAML y devuelve su contenido."""
    try:
        with open(ruta_archivo, 'r') as file:
            # Cargar el archivo YAML y devolver su contenido como un diccionario
            config = yaml.safe_load(file)
        return config
    except FileNotFoundError:
        print(f"El archivo {ruta_archivo} no se encuentra.")
    except yaml.YAMLError as e:
        print(f"Hubo un error al leer el archivo YAML: {e}")
    except Exception as e:
        print(f"Ocurrió un error inesperado: {e}")


class Timer:
    def __init__(self):
        self.start_time = None
        self.end_time = None

    def __enter__(self):
        """Se ejecuta cuando el bloque 'with' inicia."""
        self.start_time = time.time()  # Guarda el tiempo al comenzar
        return self  # Retorna la instancia para usarla dentro del 'with'

    def __exit__(self, exc_type, exc_value, traceback):
        """Se ejecuta cuando el bloque 'with' termina."""
        self.end_time = time.time()  # Guarda el tiempo al finalizar
        self.elapsed_time = self.end_time - self.start_time  # Calcula el tiempo transcurrido
        
        # Convierte el tiempo a minutos y segundos si es necesario
        minutes = self.elapsed_time // 60  # Divide por 60 para obtener los minutos
        seconds = self.elapsed_time % 60  # Obtiene el residuo para los segundos
        
        if minutes > 0:
            print(f"Tiempo de ejecución: {int(minutes)} minutos y {seconds:.4f} segundos")
        else:
            print(f"Tiempo de ejecución: {seconds:.4f} segundos")


def imprimir_mensaje(mensaje: str):
    """Imprime el mensaje en formato grande usando pyfiglet."""
    # Usamos pyfiglet para generar el arte ASCII del mensaje
    arte = pyfiglet.figlet_format(mensaje)
    print(arte)
