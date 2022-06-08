import random

from tqdm.notebook import tqdm


class AcoParada:
    # Iteraciones mínimas del algoritmo, ignorando las máximas y sin mejora
    iter_min = None
    # Iteraciones máximas del algoritmo
    iter_max = None
    # Iteraciones máximas consecutivas en las que no se mejora
    iter_sinmejoras = None
    # Ratio
    ratio_mejora = None

    def __init__(self, iter_min=None, iter_max=1000, iter_sinmejoras=None, ratio_mejora=None):
        """
        Especifica las condiciones de parada del algoritmo. Se ha de especificar al menos las iteraciones máximas o
        las iteraciones sin mejorar.
        :param iter_min: número de iteraciones mínimas garantizadas, ignorando las iteraciones máximas y sin mejorar
        :param iter_max: número de iteraciones máximas del algoritmo, ignorando los ratios de mejora
        :param iter_sinmejoras: número de iteraciones máximas consecutivas sin mejorar
        :param ratio_mejora: ratio de mejora que, en caso de ser menor, no se considerará que una solución mejore
        """
        self.iter_min = iter_min
        self.iter_max = iter_max
        self.iter_sinmejoras = iter_sinmejoras
        self.ratio_mejora = ratio_mejora

        if self.iter_max is None and self.iter_sinmejoras is None:
            raise Exception("Se ha de especificar o las iteraciones máximas totales o las iteraciones sin mejorar")
        elif self.iter_sinmejoras is not None and self.ratio_mejora is None:
            self.ratio_mejora = 0


class Aco:
    # Ajustes Fijos
    __g = None
    __parada = None, None
    __leave = True

    # Variables generadas
    __t0 = None
    __nij, __tij = None, None
    __bs, __bsc = None, None

    # Ajustes Variables
    q0 = 0.5
    alpha, beta = 1, 3
    p = 0.01
    m = 32

    def __init__(self, grafo, parada=AcoParada(), leave=True):
        """
        Inicialización del algoritmo de Colonia de Hormigas.
        :param grafo: grafo sobre el que aplicar el algoritmo
        :param parada: condiciones para la detención del algoritmo
        :param leave: si se indica como False, se ocultará la barra de progreso al acabar
        """
        self.__g = grafo
        self.__parada = parada
        self.__leave = leave

        self.__t0 = None
        self.__nij, self.__tij = None, None
        self.__bs, self.__bsc = None, None

        self.__inicializar_nij()
        self.__inicializar_tij()

    def __coste(self, sol):
        """
        Función de coste para un ACO ya instanciado.
        :param sol: solución cuyo coste se desea calcular.
        :return: el coste de la solución
        """
        return Aco.coste(self.__g, sol)

    @staticmethod
    def coste(g, sol):
        """
        Dado un grafo y una solución, calcula su coste correspondiente.
        :param g: grafo del que se extrae la solución
        :param sol: solución a la que calcular el coste
        :return: coste como decimal
        """
        s = 0
        for i in range(len(sol) - 1):
            s += g.distancia(sol[i], sol[i + 1])
        return s

    def __inicializar_nij(self):
        """
        Inicializa la matriz de costes de distancias inversa.
        """
        self.__nij = [[0.0 for _ in range(self.__g.nciudades)] for _ in range(self.__g.nciudades)]
        for i in range(self.__g.nciudades):
            for j in range(i + 1, self.__g.nciudades):
                # Dado que ir y venir tiene el mismo coste al inicio, se puede iterar como si fuese
                # una matriz triangular
                self.__nij[i][j] = 1 / self.__g.distancia(i, j)
                self.__nij[j][i] = 1 / self.__g.distancia(j, i)

    def __solucion_voraz(self):
        """
        Calcula la solución por un algoritmo voraz (búsqueda local), y devuelve la solución de este mismo
        algoritmo.
        :return: la solución obtenida por algoritmo voraz
        """
        sol = [0]
        for _ in range(1, self.__g.nciudades):
            mejor_punto, mejor_coste = None, None
            for i in range(1, self.__g.nciudades):
                if i in sol:
                    continue
                c = self.__g.distancia(sol[len(sol) - 1], i)
                if not mejor_punto or c < mejor_coste:
                    mejor_punto, mejor_coste = i, c
            sol.append(mejor_punto)
        return sol + [0]

    def __inicializar_tij(self):
        """
        Inicializa la matriz de las feromonas en base a un algoritmo de asignación voraz.
        """
        coste = self.__coste(self.__solucion_voraz())
        n = self.__g.nciudades
        self.__t0 = 1 / (coste * n)
        self.__tij = [[self.__t0 for _ in range(self.__g.nciudades)] for _ in range(self.__g.nciudades)]

    def exploracion_dirigida(self, disponibles, i):
        """
        En la asignación de caminos, esta función devolverá el siguiente nodo para una hormiga dada utilizando
        exploración dirigda. Esta consiste en asignar unas proabilidades, y en base a una distribución ordenada,
        decidir de manera aleatoria por probabilidad a cual se avanza.
        :param disponibles: ciudades pendientes por visitar
        :param i: ciudad actual
        :return: ciudad a visitar a continuación
        """

        # Se calcula la suma de todas las probabilidades
        s = sum((self.__tij[i][l] ** self.alpha) * (self.__nij[i][l] ** self.beta) for l in disponibles)
        opciones = []
        for l in disponibles:
            # Para cada opción, calcular su probabilidad divida por la total
            opciones.append((l, ((self.__tij[i][l] ** self.alpha) * (self.__nij[i][l] ** self.beta)) / s))
        # Se ordena las soluciones por probabilidad de mayor a menor
        opciones = sorted(opciones, key=lambda o: 1 / o[1])

        # Se genera un número aleatorio de "corte"
        m = random.uniform(0, 1)
        j, agregado = None, 0
        for opcion in opciones:
            # Se suma la probabilidad agregada, y se almacena la ciudad a visitar
            agregado += opcion[1]
            j = opcion[0]
            # Si la suma agregada es mayor que el corte, se para de iterar y se devuelve la ciudad encontrada
            if agregado >= m:
                break
        return j

    def explotacion(self, disponibles, i):
        """
        En la asignación de caminos, la explotación consiste en realizar una búsqueda local para encontrar
        el mejor camino. Se tiene en cuenta la distancia y la feromona, pero sólo se tendrá en cuenta cual es el
        próximo mejor en base a esos dos parámetros.
        :param disponibles: ciudades pendientes de visitar
        :param i: ciudad actual
        :return: ciudad que se ha de visitar
        """
        opciones = []
        for l in disponibles:
            # Para cada ciudad disponible, almacenar un valor indicando feromona y distancia
            opciones.append((l, (self.__tij[i][l] ** self.alpha) * (self.__nij[i][l] ** self.beta)))
        opciones = sorted(opciones, key=lambda o: o[1])
        # Se devuelve la ciudad cuyo valor sea el más alto
        j = opciones[len(opciones) - 1][0]
        return j

    def calcular_soluciones_k(self):
        """
        Para cada hormiga de la colonia, se calcula una solución válida. Para esto, hay dos opciones de calcular los
        nodos siguientes: mediante exploración dirigida o explotación. Esto se hace de manera aleatoria, con una
        probabilidad dada por el atributo q0.
        :return: lista de soluciones, una por hormiga
        """
        sols = []
        for k in range(self.m):
            sol = [0]
            # Para cada hormiga, se almacena en una lista las ciudades pendientes de visitar
            disponibles = list(range(1, self.__g.nciudades))
            # Se itera mientras queden ciudades por visitar
            while len(disponibles) > 0:
                # La ciudad i será la última visitada
                i = sol[len(sol) - 1]
                q = random.uniform(0, 1)
                # En función de un valor aleatorio, se decide que la ciudad a visitar j será extraida por
                # explotación o por exploración dirigida
                if q < self.q0:
                    j = self.explotacion(disponibles, i)
                else:
                    j = self.exploracion_dirigida(disponibles, i)
                # Se inserta el elemento en la solución de la hormiga y se elimina de las disponibles
                sol.append(j)
                disponibles.remove(j)
            # Se guarda la solución en la lista de soluciones, una por hormiga
            sols.append(sol + [0])
        return sols

    def seleccionar_mejor_solucion(self, sols):
        """
        Dada una lista de soluciones obtenidas por el algoritmo, se extrae la mejor en función de su coste. Se
        tiene en cuenta también la mejor solución global, por lo que es posible que se devuelva la misma solución
        global.
        :param sols: lista de soluciones
        :return: mejor solución, incluyendo en el cálculo la global
        """
        mejor_solucion, mejor_coste = None, None
        for sol in sols:
            coste = self.__coste(sol)
            # Si se mejora el coste con respecto a otras soluciones locales, guardar como mejor solución local
            if mejor_solucion is None or coste < mejor_coste:
                mejor_solucion, mejor_coste = sol, coste

        # Si se mejora el coste de la mejor solución global, devolver la nueva
        if self.__bs is None or mejor_coste < self.__bsc:
            return mejor_solucion, mejor_coste
        # Sino, devolver la global
        return self.__bs, self.__bsc

    def actualizar_feromona_global(self):
        """
        Cogiendo la mejor solución global, se aplica la fórmula de actualización de feromonas a los caminos
        visitados por esta solución (sólo de ida, no en ida y vuelta ya que los caminos tienen solución). El
        resultado se almacena en la matriz de feromonas.
        """
        for k in range(1, len(self.__bs)):
            i, j = self.__bs[k - 1], self.__bs[k]
            self.__tij[i][j] = (1 - self.p) * self.__tij[i][j] + self.p * (1 / self.__bsc)

    def solucion(self):
        """
        Se devuelve la solución del camino óptimo de un grafo mediante el algoritmo de colonias de hormigas.
        Dependiendo de la condición de parada, se llamará a una estructura del algoritmo u otra, siendo la por
        defecto la de un máximo de 1000 iteraciones. Además, al acabar de encontrar una solución, se resetean todas
        las matrices temporales, permitiendo ejecutarlo de nuevo cambiando los atributos públicos.
        :return: solución óptima encontrada
        """
        iters = 0
        iters_sinmejora = 0

        # Iniciar la abrra de progreso en función de la condición máxima de parada
        iters_max = None
        if self.__parada.iter_max is not None:
            iters_max = self.__parada.iter_max
        elif self.__parada.iter_sinmejoras is not None:
            iters_max = self.__parada.iter_sinmejoras
        pbar = tqdm(total=iters_max, leave=self.__leave)

        while True:
            if self.__parada.iter_max is not None:
                iters += 1
                pbar.update()

            # Se calculan las soluciones para cada hormiga
            sols = self.calcular_soluciones_k()

            # Se selecciona la mejor solución de las encontradas, incluyendo la actual
            bs, bsc = self.seleccionar_mejor_solucion(sols)
            # Y se calcula el ratio de mejora con respecto a la anterior, actualizando la global tras esto
            if self.__bs is None:
                ratio_mejora = 1
            else:
                ratio_mejora = 1 - (bsc / self.__bsc)
            self.__bs, self.__bsc = bs, bsc

            # Se actualiza la feromona en base a esta solución encontrada
            self.actualizar_feromona_global()

            # Si se especifica una condición de parada de iteraciones sin mejora, resetear el contador si el ratio
            # de mejora es superior al especificado (por defecto, es 0)
            if self.__parada.iter_sinmejoras is not None and ratio_mejora > self.__parada.ratio_mejora:
                iters_sinmejora = 0
            # y en caso contrario, incrementar el contador
            elif self.__parada.iter_sinmejoras is not None:
                iters_sinmejora += 1
            # Actualizar la barra de progreso si la condición de parada máxima es la sin mejoras
            if self.__parada.iter_max is None:
                pbar.update()

            # Si se indica una conidición de iteraciones mínimas, seguir ejecutando el algoritmo pase lo que pase
            if self.__parada.iter_min is not None and iters < self.__parada.iter_min:
                continue
            # Si se especifica una condición de parada con ciertas iteraciones sin mejorar y se alcanzan, detener
            if self.__parada.iter_sinmejoras is not None and iters_sinmejora >= self.__parada.iter_sinmejoras:
                break
            # Si se especifica una condición de parada con un máximo de iteraciones y se alcanzan, detener
            if self.__parada.iter_max is not None and iters >= self.__parada.iter_max:
                break
        # Cerrar la barra de progreso, forzando el completado
        if self.__parada.iter_sinmejoras is not None and pbar.n < pbar.total:
            pbar.update(pbar.total - pbar.n)
        pbar.close()

        # Almacenar la solución actual en una tupla
        sol = self.__bs, self.__bsc

        # Resetear las matrices de cálculo y la solución obtenida
        self.__init__(self.__g, self.__parada, self.__leave)

        # Devolver la solución
        return sol
