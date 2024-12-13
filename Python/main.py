from src.helpers import load_config, Timer, imprimir_mensaje


if __name__ == "__main__":

    with Timer() as t:

        # Cargamos el archivo config.yaml
        config = load_config("config.yaml")
        

        



        # Imprimimos el mensaje de finalización de la ejecución
        imprimir_mensaje("Todo bien !")
