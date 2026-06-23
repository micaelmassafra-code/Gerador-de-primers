

import requests

# --------------------------------------------
# 1. Cálculo da porcentagem de GC
# --------------------------------------------
def calc_gc(seq):
    seq = seq.upper()
    if len(seq) == 0:
        return 0
    gc = sum(1 for base in seq if base in "GC")
    return (gc / len(seq)) * 100


# --------------------------------------------
# 2. Cálculo da Temperatura de Melting (Tm)
# --------------------------------------------
def calc_tm(seq):
    seq = seq.upper()
    at = sum(1 for base in seq if base in "AT")
    gc = sum(1 for base in seq if base in "GC")
    return 2 * at + 4 * gc


# --------------------------------------------
# 3. Verificação do Clamp GC
# --------------------------------------------
def has_gc_clamp(seq):
    seq = seq.upper()
    return seq.endswith("G") or seq.endswith("C")


# --------------------------------------------
# 4. Checagem de homopolímeros
# --------------------------------------------
def check_homopolymer(seq, limit=3):
    seq = seq.upper()
    bases = ["A", "T", "C", "G"]
    for b in bases:
        if b * (limit + 1) in seq:
            return False
    return True


# --------------------------------------------
# 5. Reverse Complement (manual)
# --------------------------------------------
def reverse_complement(seq):
    seq = seq.upper()
    complemento = {"A": "T", "T": "A", "C": "G", "G": "C"}
    return "".join(complemento[b] for b in reversed(seq))


# --------------------------------------------
# 6. Verificação de primer válido
# --------------------------------------------
def primer_valido(seq, gc_min=40, gc_max=60, tm_min=55, tm_max=65, homo_limit=3):
    gc = calc_gc(seq)
    tm = calc_tm(seq)
    clamp = has_gc_clamp(seq)
    homo = check_homopolymer(seq, limit=homo_limit)

    return (gc_min <= gc <= gc_max) and (tm_min <= tm <= tm_max) and clamp and homo


# --------------------------------------------
# 7. Geração de primers válidos a partir da sequência
# --------------------------------------------
def gerar_primers(seq, min_size=18, max_size=24):
    seq = seq.upper().replace(" ", "").replace("\n", "")
    primers = []

    for i in range(len(seq)):
        for size in range(min_size, max_size + 1):
            cand = seq[i:i+size]
            if len(cand) < size:
                continue

            if primer_valido(cand):
                primers.append({
                    "seq": cand,
                    "gc": calc_gc(cand),
                    "tm": calc_tm(cand),
                    "pos": i
                })

    return primers


# --------------------------------------------
# 8. Verificação de pares Forward/Reverse
# --------------------------------------------
def verificar_pares(seq, primers_f, primers_r, min_amp=100, max_amp=200, max_delta_tm=5):
    seq = seq.upper()
    pares = []

    for p1 in primers_f:
        for p2 in primers_r:
            pos_fwd = p1["pos"]
            rev_comp = reverse_complement(p2["seq"])
            pos_rev = seq.find(rev_comp)

            if pos_rev == -1 or pos_fwd >= pos_rev:
                continue

            amplicon = pos_rev - pos_fwd
            delta_tm = abs(p1["tm"] - p2["tm"])

            if (min_amp <= amplicon <= max_amp) and (delta_tm <= max_delta_tm):
                pares.append({
                    "forward": p1["seq"],
                    "reverse": p2["seq"],
                    "amplicon": amplicon,
                    "tm_fwd": p1["tm"],
                    "tm_rev": p2["tm"],
                    "delta_tm": delta_tm
                })

    return pares


# ------------------------------------------------------------
# 9. Gerar pares forçados (corrigido)
# ------------------------------------------------------------
def pares_forcados(fwds, revs):
    forcados = []
    top_f = fwds[:10]
    top_r = revs[:10]

    for pf in top_f:
        for pr in top_r:
            delta_tm = abs(pf["tm"] - pr["tm"])
            forcados.append({
                "Forward": pf["seq"],
                "Reverse": pr["seq"],
                "Tm_F": pf["tm"],
                "Tm_R": pr["tm"],
                "Delta_Tm": delta_tm
            })

    forcados.sort(key=lambda x: x["Delta_Tm"])
    return forcados[:5]


# ------------------------------------------------------------
# 10. Buscar primer individual no NCBI (corrigido)
# ------------------------------------------------------------
def buscar_ncbi(primer):
    primer = primer.upper().strip()

    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
    params = {
        "db": "nucleotide",
        "term": f"\"{primer}\"[Sequence]",
        "retmode": "json"
    }

    resp = requests.get(url, params=params)

    if resp.status_code != 200:
        print("Erro ao acessar NCBI:", resp.status_code)
        return

    dados = resp.json()
    total = dados["esearchresult"]["count"]

    print("\n===== RESULTADO NCBI =====")
    print(f"Primer buscado: {primer}")
    print(f"Quantidade de sequências no GenBank que CONTÊM esse primer: {total}")

    if total != "0":
        print("IDs encontrados:", dados["esearchresult"]["idlist"][:10])
    else:
        print("Nenhuma sequência contém este primer.")


# ------------------------------------------------------------
# Pipeline completo
# ------------------------------------------------------------
def executar_pipeline(seq):
    print("\nGerando primers...")
    primers_f = gerar_primers(seq)
    primers_r = gerar_primers(reverse_complement(seq))

    print("\n--- 10 MELHORES PRIMERS FORWARD ---")
    for i, p in enumerate(primers_f[:10], 1):
        print(f"{i}. {p['seq']} | GC={p['gc']:.2f}% | Tm={p['tm']}")

    print("\n--- 10 MELHORES PRIMERS REVERSE ---")
    for i, p in enumerate(primers_r[:10], 1):
        print(f"{i}. {p['seq']} | GC={p['gc']:.2f}% | Tm={p['tm']}")

    pares = verificar_pares(seq, primers_f, primers_r)
    print("\n--- PARES COMPATÍVEIS (100–200 pb) ---")

    if pares:
        for p in pares[:5]:
            print("-" * 40)
            print(p)
    else:
        print("Nenhum par compatível encontrado!")

    input("\nPressione ENTER para voltar ao menu...")


# ------------------------------------------------------------
# MENU
# ------------------------------------------------------------
def menu():
    print("\n===== MENU =====")
    print("1 - Calcular Tm")
    print("2 - Calcular GC%")
    print("3 - Complemento reverso")
    print("4 - Rodar pipeline completo")
    print("5 - Gerar pares forçados")
    print("6 - Buscar primer de pares compatíveis no NCBI")
    print("7 - Buscar primer de pares forçados no NCBI")
    print("0 - Sair")
    return input("Escolha: ")


# ------------------------------------------------------------
# PROGRAMA PRINCIPAL
# ------------------------------------------------------------
print("===============================================")
print("    SISTEMA DE ANÁLISE E SELEÇÃO DE PRIMERS    ")
print("                Projeto IS6110                  ")
print("===============================================")

seq_input = input("\nDigite a sequência de DNA a ser analisada:\n→ ").upper()
linhas = [l for l in seq_input.splitlines() if not l.startswith(">")]
seq = "".join(linhas).replace(" ", "").replace("\t", "").upper()

pares_forcados_salvos = []

while True:
    op = menu()

    if op == "1":
        p = input("Primer: ").upper()
        print("Tm =", calc_tm(p))
        input("\nPressione ENTER para voltar ao menu...")

    elif op == "2":
        p = input("Primer: ").upper()
        print("GC% =", calc_gc(p))
        input("\nPressione ENTER para voltar ao menu...")

    elif op == "3":
        p = input("Sequência: ").upper()
        print("RC =", reverse_complement(p))
        input("\nPressione ENTER para voltar ao menu...")

    elif op == "4":
        executar_pipeline(seq)

    elif op == "5":
        print("\nGerando primers forçados...")
        fwd = gerar_primers(seq)
        rev = gerar_primers(reverse_complement(seq))
        pares_forcados_salvos = pares_forcados(fwd, rev)
        print("\n--- 5 MELHORES PARES FORÇADOS ---\n")
        for p in pares_forcados_salvos:
            print(p)
        input("\nPressione ENTER para voltar ao menu...")

    elif op == "6":
        primers_f = gerar_primers(seq)
        primers_r = gerar_primers(reverse_complement(seq))
        pares = verificar_pares(seq, primers_f, primers_r)

        if not pares:
            print("Nenhum par compatível encontrado!")
        else:
            for par in pares[:5]:
                print("\nPar selecionado:")
                print(par)

                print("\nBuscar qual primer no NCBI?")
                print("1 - Forward")
                print("2 - Reverse")
                escolha = input("→ ")

                if escolha == "1":
                    buscar_ncbi(par["forward"])
                elif escolha == "2":
                    buscar_ncbi(par["reverse"])

        input("\nPressione ENTER para voltar ao menu...")

    elif op == "7":
        if not pares_forcados_salvos:
            print("Nenhum par forçado gerado ainda.")
            input("\nPressione ENTER para voltar ao menu...")
            continue

        for par in pares_forcados_salvos:
            print("\nPar selecionado:")
            print(par)

            print("\nBuscar qual primer no NCBI?")
            print("1 - Forward")
            print("2 - Reverse")
            escolha = input("→ ")

            if escolha == "1":
                buscar_ncbi(par["Forward"])
            elif escolha == "2":
                buscar_ncbi(par["Reverse"])

        input("\nPressione ENTER para voltar ao menu...")

    elif op == "0":
        print("Encerrando.")
        break

    else:
        print("Opção inválida!")
        input("\nPressione ENTER para voltar ao menu...")
