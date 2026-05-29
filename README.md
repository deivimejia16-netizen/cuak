# cuak
nada interesante
import matplotlib.pyplot as plt

def calcular_porcentajes_individuales(cadena_adn):
    """
    Calcula los porcentajes globales de cada una de las 4 bases individuales (A, T, G, C).
    """
    longitud = len(cadena_adn)
    if longitud == 0:
        return 0, 0, 0, 0
    
    p_a = (cadena_adn.count('A') / longitud) * 100
    p_t = (cadena_adn.count('T') / longitud) * 100
    p_g = (cadena_adn.count('G') / longitud) * 100
    p_c = (cadena_adn.count('C') / longitud) * 100
    
    return p_a, p_t, p_g, p_c

def analizar_ventanas_gc(secuencia, tamano_ventana=100):
    """
    Analiza la secuencia en ventanas fijas de 100 nt calculando únicamente el %GC local.
    """
    posiciones = []
    valores_gc_local = []
    
    for i in range(0, len(secuencia), tamano_ventana):
        ventana = secuencia[i:i + tamano_ventana]
        
        # para que lea marcos de 100 nt
        if len(ventana) == tamano_ventana:
            longitud_v = len(ventana)
            g = ventana.count('G')
            c = ventana.count('C')
            gc_local = ((g + c) / longitud_v) * 100
            
            valores_gc_local.append(gc_local)
            posiciones.append(i)
            
    return posiciones, valores_gc_local

def graficar_perfil_gc(posiciones, gc_locales, gc_global):
    """
    Genera una gráfica enfocada únicamente en el perfil de %GC (Guanina + Citosina).
    """
    plt.figure(figsize=(12, 6))
    
    # Grafica de %GC
    plt.plot(posiciones, gc_locales, color='#2c3e50', linewidth=1.8, label='%GC Local (Guanina + Citosina)')
    
    # pomedio de y punto de referencia %GC Global
    plt.axhline(y=gc_global, color='#e74c3c', linestyle='--', linewidth=1.5, 
                label=f'Promedio Global %GC ({gc_global:.2f}%)')
    
    # grafica variables organizacion
    plt.title('Perfil de Composición Genómica (%GC) - Arquea Extremófila', fontsize=13, fontweight='bold')
    plt.xlabel('Posición en la Secuencia de ADN (pb)', fontsize=11)
    plt.ylabel('Porcentaje (%) de GC', fontsize=11)
    plt.ylim(-5, 105)
    plt.grid(True, linestyle=':', alpha=0.5)
    plt.legend(loc='best', fontsize=10)
    
    # muestra la grafica en imagen
    plt.savefig('perfil_gc_unico.png', dpi=300, bbox_inches='tight')
    print("\n[INFO] Gráfica de %GC exportada exitosamente como 'perfil_gc_unico.png'.")
    plt.show()

# chambeo
if __name__ == "__main__":
    print("=== SCRIPT DE ANÁLISIS GENÓMICO INTERACTIVO ===")
    
    # secuencia de entrada
    print("\nPor favor, introduce o pega la secuencia de ADN (A, T, G, C) y presiona ENTER:")
    entrada_usuario = input()
    
    # flitro de ruido
    secuencia_adn = "".join([char for char in entrada_usuario if char.upper() in 'ATGC']).upper()
    longitud_total = len(secuencia_adn)
    
    if longitud_total == 0:
        print("[ERROR] No se detectó ninguna secuencia nucleotídica válida.")
    else:
        print(f"\n[ÉXITO] Secuencia procesada correctamente.")
        print(f"Longitud útil total: {longitud_total} pares de bases (pb).")
        
        # 1. Cálculo de métricas globales individuales (Desglose de bases solicitado)
        pct_a, pct_t, pct_g, pct_c = calcular_porcentajes_individuales(secuencia_adn)
        gc_global = pct_g + pct_c
        at_global = pct_a + pct_t
        
        print("\n" + "="*45)
        print("          MÉTRICAS GENÓMICAS GLOBALES        ")
        print("="*45)
        print(f" Adenina (A) global : {pct_a:.2f}%")
        print(f" Timina (T) global  : {pct_t:.2f}%")
        print(f" Guanina (G) global : {pct_g:.2f}%")
        print(f" Citosina (C) global: {pct_c:.2f}%")
        print("-"*45)
        print(f" TOTAL CONTENIDO %GC: {gc_global:.2f}%")
        print(f" TOTAL CONTENIDO %AT: {at_global:.2f}%")
        print("="*45)
        
        # 2. Análisis local por ventanas de 100 nt
        posiciones, gc_locales = analizar_ventanas_gc(secuencia_adn, tamano_ventana=100)
        
        # 2.1 Mostrar las densidades locales calculadas en la terminal
        print("\n--- Densidad de %GC por Ventanas (Primeras posiciones) ---")
        max_vistas = min(10, len(posiciones))
        for i in range(max_vistas):
            pos = posiciones[i]
            print(f"Ventana Rango: pb {pos:04d} - {pos+100:04d} | %GC Local: {gc_locales[i]:.2f}%")
            
        # 3. Grafica de %GC
        if len(posiciones) > 0:
            graficar_perfil_gc(posiciones, gc_locales, gc_global)
        else:
            print("\n[AVISO] La secuencia es menor a 100 pb; no se pueden generar ventanas fijas para graficar.")
