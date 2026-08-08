# Detector jerárquico Wigner-star multiqubit

Versión española estilo arXiv preparada para revisión de Enrique E. Casanova
Benítez y Oscar Jasel Berra Montiel.

Repositorio público asociado al manuscrito:
<https://github.com/enrik03/wigner-star-multiqubit-entanglement>

## Artefactos

- `main.tex`: manuscrito fuente.
- `main.pdf`: PDF compilado de 12 paginas.
- `references.bib`: 27 referencias citadas, sin entradas huerfanas.
- `scripts/verify_detector_jerarquico.py`: cinco casos de tres qubits, umbrales y SDP.
- `scripts/verify_four_qubit_extension.py`: GHZ4 en siete biparticiones.
- `scripts/validate_results.py`: puerta numerica independiente.
- `scripts/check_manuscript_sync.py`: sincronizacion entre JSON, tablas, citas y TeX.
- `results/*.json`: salidas reproducibles y puertas de integridad.
- `AUDITORIA_PARA_JASEL.md`: evaluacion critica y trabajo pendiente.
- `CITATION.cff`: metadatos para citar este paquete reproducible.

## Entorno

```powershell
python -m pip install -r requirements.txt
```

La implementacion usa NumPy y CVXPY. CLARABEL es el solver primario y SCS queda
como respaldo.

## Verificación

```powershell
python -m py_compile scripts/verify_detector_jerarquico.py `
  scripts/validate_results.py `
  scripts/verify_four_qubit_extension.py `
  scripts/check_manuscript_sync.py

python scripts/verify_detector_jerarquico.py
python scripts/validate_results.py
python scripts/verify_four_qubit_extension.py
python scripts/check_manuscript_sync.py
```

Resultados de la entrega actual:

- puerta de tres qubits: `passed=true`;
- error maximo Wigner-star frente a matrices: `9.094947017729282e-13`;
- puerta GHZ4: `passed=true`, siete cortes, error de momentos cero;
- sincronizacion manuscrito-resultados: `passed=true`;
- compilacion Tectonic: exitosa, sin avisos LaTeX, citas indefinidas ni cajas desbordadas.

## Cita

GitHub puede generar una cita desde `CITATION.cff`. Para una versión sometida o
publicada, debe citarse además el artículo y el DOI archivado correspondiente.

## Compilacion

Con Tectonic:

```powershell
tectonic -X compile main.tex --keep-logs --keep-intermediates
```

El aviso externo de Fontconfig observado en Windows pertenece al runtime y no
al documento; `main.log` no contiene avisos de LaTeX.
