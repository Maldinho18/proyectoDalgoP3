import random
import os

def generate_test_case(n=10**4, k=10**3, alphabet='ACGT', out_filename="large_test_case.txt"):
    """
    Genera un test con:
      - 1 caso de prueba
      - n k-mers de longitud k
      - Salida en el fichero out_filename en el directorio actual.
    """
    # Crear superstring S de longitud L = n + k - 1
    L = n + k - 1
    S = ''.join(random.choices(alphabet, k=L))
    # Extraer todos los k-mers solapados
    kmers = [S[i:i+k] for i in range(L - k + 1)]
    random.shuffle(kmers)

    # Escribir input
    with open(out_filename, "w") as f:
        f.write("1\n")
        f.write(f"{len(kmers)} {k}\n")
        for mer in kmers:
            f.write(mer + "\n")
    # Escribir expected output por separado
    exp_fn = "expected_output.txt"
    with open(exp_fn, "w") as f:
        f.write(S + "\n")

    print(f"Test generado: n={len(kmers)}, k={k}")
    print(f"- Input:  {out_filename}")
    print(f"- Output: {exp_fn}")

if __name__ == "__main__":
    # Fija semilla para reproducibilidad
    random.seed(42)
    generate_test_case()
