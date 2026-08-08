# Auditoria para Jasel

**Manuscrito:** *Detector jerarquico de entrelazamiento multiqubit desde
funciones de Wigner Stratonovich-Weyl*
**Fecha de corte:** 8 de agosto de 2026
**Version auditada:** `main.tex` / `main.pdf` de esta carpeta

## Veredicto

**Puntuacion interna como preprint cientifico para mostrar a Jasel: 95/100.**

**Preparacion inmediata para enviar a EPJ Plus: 91/100.** La diferencia no es
matematica: falta aprobacion del coautor, version inglesa/formato editorial,
repositorio archivado con DOI y una ultima lectura humana especializada.

Esta es una auditoria interna multienfo (editorial, metodologia, dominio y
abogado del diablo); no sustituye un arbitraje externo.

## Desglose de la puntuacion cientifica

| Dimension | Puntos | Evidencia |
|---|---:|---|
| Correccion matematica | 20/20 | Normalizacion SW fija; regla de traza; Newton completo; prueba de positividad del mapa; SDP equivalente. |
| Novedad y posicionamiento | 17/20 | Sintesis Wigner-star de tres capas y uso de `lambda_min(H_l)`; los criterios base son conocidos. |
| Ejemplos y controles | 14/15 | GHZ3, W3, GHZ4, dos familias ruidosas y mezcla biseparable NPT en todos los cortes. |
| Reproducibilidad | 15/15 | Tres verificadores, JSON, puerta independiente y sincronizacion automatica con el TeX. |
| Cobertura bibliografica | 10/10 | 27/27 referencias citadas, incluidas fuentes primarias de 2025-2026. |
| Claridad y limites | 10/10 | Se separan NPT por corte, certificado GME suficiente y exclusion de mezclas PPT. |
| Ajuste editorial | 9/10 | PDF limpio de 12 paginas; falta convertir la entrega espanola al formato final de revista. |

## Contribucion defendible

1. El dato primario es el simbolo de Wigner SW y toda la cadena finita se
   escribe como `W -> reflexion -> potencias star -> momentos -> espectro`.
2. `F_star^(S)` es exacto para PPT/NPT de cada corte cuando se usan los `2^N`
   momentos necesarios. No se presenta como monotono LOCC ni medida universal.
3. La segunda capa traslada mapas GME a fase y usa la condicion completa
   `H_l >= 0`. El funcional basado en `-lambda_min(H_l)` es mas fuerte que mirar
   solo `det(H_l)` y evita fallos de paridad o rango.
4. La tercera capa formula en el cono Wigner-star el SDP de testigos
   completamente descomponibles. `N_GME^star > 0` equivale a estar fuera de las
   mezclas PPT y, por tanto, certifica GME.
5. El contraejemplo biseparable Bell demuestra que NPT en todos los cortes no
   basta para GME mixta. El caso GHZ4 verifica la extension sobre siete cortes,
   no solo sobre los tres cortes tripartitos.

## Mejoras realizadas en esta version

- Anclaje explicito en el producto estrella de qubits de Berra-Montiel,
  Molgado y Sanchez-Cordova.
- Normalizacion general de la traza del mapa:
  `T_N[M_T^(N) W] = |B_N| + 2^N c_N`; da 11 para N=3 y 55 para N=4.
- Lema autosuficiente que demuestra por que
  `c_N=(2^(N-1)-2)/2` hace positivo el mapa sobre todo el casco biseparable.
- GHZ4 introducido directamente por sus coeficientes Wigner; siete cortes con
  `F_star=-1/4`, mapa minimo `-1/2`, primera deteccion en `H_2` y
  `N_GME^star=0.499999993`.
- Referencia directa a Zhang et al. sobre momentos PT y su extension GHZ/W.
- Reclamo Ha-Kye corregido: existen estados GME dentro de las mezclas PPT; no se
  afirma PPT simultanea en todos los cortes sin evidencia.
- Puerta `check_manuscript_sync.py` para impedir divergencias entre tablas,
  umbrales, citas y salidas JSON.

## Evidencia reproducible

- Tres qubits: `results/detector_jerarquico_checks.json`.
- GHZ4: `results/four_qubit_extension.json`.
- Integridad numerica: `results/manuscript_sync_check.json` y salida de
  `scripts/validate_results.py`.
- Error maximo entre momentos estrella y trazas matriciales independientes:
  `9.094947017729282e-13`.
- GHZ4: error de coeficientes, mapa y momentos igual a cero en doble precision.
- PDF: 12 paginas, 27 referencias, autores y acentos extraibles, sin caracteres
  de reemplazo, referencias indefinidas ni desbordamientos.

## Debilidades que deben permanecer visibles

1. **Novedad incremental.** La contribucion principal es una sintesis exacta en
   fase, la condicion espectral completa de Hankel y una logica jerarquica. No
   debe venderse como un criterio de separabilidad universal nuevo.
2. **PPT-GME.** Los mapas de transposicion son descomponibles y el SDP de
   mezclas PPT no detecta todos los estados GME. Cero siempre significa
   inconcluso, salvo en clases especiales demostradas.
3. **Dependencia del mapa.** La eleccion estandar/modificada usa conocimiento de
   familias GHZ/W. Una optimizacion sobre mapas admisibles sigue abierta.
4. **Costo experimental.** La tabla completa de W tiene `4^N-1` parametros. El
   manuscrito no demuestra una ventaja de muestras frente a tomografia; solo
   identifica donde momentos multicopia podrian evitarla.
5. **Ruido estadistico.** Las pruebas actuales usan aritmetica numerica ideal.
   Faltan intervalos certificados para momentos experimentales cerca de la
   frontera PPT.
6. **Solver GHZ4.** CLARABEL reporta `optimal_inaccurate` para el SDP de 16x16,
   aunque el valor difiere de 1/2 por menos de `7e-9` y las capas espectral y
   Hankel proporcionan certificados independientes.

## Lo que Jasel debe decidir

1. Aprobar que la novedad se formule como **arquitectura Wigner-star
   jerarquica**, no como nueva medida universal.
2. Confirmar la convencion SW y la compatibilidad exacta con su preprint de
   producto estrella de 2026.
3. Decidir si el envio inicial sera esta version arXiv espanola o una version
   inglesa simultanea.
4. Aprobar autoria, contribuciones, afiliaciones y declaraciones.
5. Archivar una version estable del repositorio publico con DOI antes del envio.

## Camino corto a EPJ Plus

1. Congelar primero esta version tecnica con Jasel.
2. Traducir al ingles cientifico sin cambiar ecuaciones ni reclamos.
3. Pasar al formato Springer/EPJ Plus y preparar carta de presentacion.
4. Ejecutar un segundo solver de alta precision para GHZ4 o reportar residuos
   primal/dual en un suplemento.
5. Archivar codigo/JSON con DOI y repetir las cuatro puertas desde un entorno
   limpio.

## Fuentes primarias revisadas en linea

- [Berra-Montiel, Molgado y Sanchez-Cordova, arXiv:2604.05170](https://arxiv.org/abs/2604.05170)
- [Mukherjee et al., arXiv:2506.00162 / PRA 112, 062428](https://arxiv.org/abs/2506.00162)
- [Jungnitsch, Moroder y Guhne, arXiv:1010.6049](https://arxiv.org/abs/1010.6049)
- [Ha y Kye, arXiv:1512.04693 / PRA 93, 032315](https://arxiv.org/abs/1512.04693)
- [Zhang et al., arXiv:2404.19308 / Annalen der Physik 534, 2200289](https://arxiv.org/abs/2404.19308)

## Recomendacion final

Esta version ya es adecuada para una reunion tecnica con Jasel y para circular
como preprint interno. No recomiendo ampliar mas ejemplos antes de su lectura:
el siguiente incremento de valor no viene de otra tabla GHZ/W, sino de elegir
con el coautor entre dos rutas: un mapa indecomponible para cubrir parte de
PPT-GME o un protocolo experimental con costo y varianza cuantificados.
