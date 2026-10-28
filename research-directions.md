# Research Directions

Seven projects that connect the architecture's three ideas to existing research lines, ordered roughly from easiest to start to hardest. The oracle memory is what the literature calls quantum random access memory (QRAM). The classical handling of X and SWAP is Pauli frame tracking and qubit relabeling. The active/inactive split matches zoned hardware.

1. **Build the gate-absorption layer and measure it.** Write the rules for how an X flag and an address swap propagate through every other gate. Flags pass cleanly through Clifford gates but not through T gates, which is the real limit of the idea. Compile the Fourier transform, Grover and Shor circuits in Qiskit and count the gates removed. Start from Knill 2005, Riesebos 2017, McKay 2017 and the SABRE routing paper.
2. **Define the classical controller's instruction set.** Specify instructions for remapping addresses, setting flags, and firing operations on a tick, then write a small simulator. The eQASM and QuMA papers show how a working controller does timing. Pure software; can start now.
3. **Cost the oracle lookup.** The paper treats a memory lookup as free, and this is its largest gap. Work out gate count and depth per lookup for the bucket-brigade design, table-lookup circuits, and the select-swap variant. Redo the Grover qubit table with lookup cost included. Read Jaques and Rattew first.
4. **Test the shared-oracle claim against parallel queries.** A standard bucket-brigade memory answers one query at a time, so three Grover instances would have to take turns. Fat-Tree QRAM addresses exactly this. Zalka's result also matters: parallel Grover instances give only a square-root speedup. Deliverable: a curve of qubits against query throughput.
5. **Replace the hypothetical error operator with a real noise model.** Use amplitude damping and dephasing channels (Nielsen & Chuang ch. 8), with idle errors on the inactive lines. Simulate in Stim or Qiskit Aer. Key question: the best number of operations per tick, since packing more in raises crosstalk.
6. **Map the architecture onto zoned hardware.** Trapped-ion machines separate memory zones from interaction zones; neutral-atom machines have storage, entangling and readout zones. In both, a swap by relabeling corresponds to physically moving a qubit. Thaker et al. 2006 is the closest prior architecture.
7. **Run a full resource estimate.** Check the scalability objective with physical qubit counts and runtimes, with and without the architecture. Follow Gidney and Ekerå, or use the Azure Quantum Resource Estimator (Beverland et al.).

Ideas 1 and 2 need only a laptop. Ideas 3 and 4 decide whether the central claim survives. The ISEF research plan scopes the project to ideas 1, 3, 4 and 5.

## Papers (confirmed against arXiv, 2026-10-01)

| Paper | arXiv | Idea |
|---|---|---|
| Giovannetti, Lloyd, Maccone, "Quantum random access memory", 2008 | 0708.1879 | 3 |
| Giovannetti, Lloyd, Maccone, "Architectures for a quantum random access memory", 2008 | 0807.4994 | 3 |
| Babbush et al., "Encoding Electronic Spectra in Quantum Circuits with Linear T Complexity", 2018 | 1805.03662 | 3 |
| Low, Kliuchnikov, Schaeffer, "Trading T gates for dirty qubits", 2018 | 1812.00954 | 3 |
| Jaques, Rattew, "QRAM: A Survey and Critique", 2023 | 2305.10310 | 3, 4 |
| Di Matteo, Gheorghiu, Mosca, "Fault tolerant resource estimation of quantum random-access memories", 2019 | 1902.01329 | 3 |
| Arunachalam et al., "On the robustness of bucket brigade quantum RAM", 2015 | 1502.03450 | 3, 5 |
| Hann, Lee, Girvin, Jiang, "Resilience of quantum random access memory to generic noise", 2020 | 2012.05340 | 5 |
| Xu, Hann, Foxman, Girvin et al., "Systems Architecture for Quantum Random Access Memory", 2023 | 2306.03242 | 4 |
| Xu, Lu, Ding, "Fat-Tree QRAM", 2025 | 2502.06767 | 4 |
| Zalka, "Grover's quantum searching algorithm is optimal", 1997 | quant-ph/9711070 | 4 |
| McKay et al., "Efficient Z-Gates for Quantum Computing", 2016 | 1612.00858 | 1 |
| Li, Ding, Xie, "Tackling the Qubit Mapping Problem for NISQ-Era Quantum Devices", 2018 | 1809.02573 | 1 |
| Fu et al., "An Experimental Microarchitecture for a Superconducting Quantum Processor", 2017 | 1708.07677 | 2 |
| Fu et al., "eQASM: An Executable Quantum Instruction Set Architecture", 2018 | 1808.02449 | 2 |
| Fowler et al., "Surface codes: Towards practical large-scale quantum computation", 2012 | 1208.0928 | 5 |
| Terhal, "Quantum Error Correction for Quantum Memories", 2013 | 1302.3428 | 5 |
| Google Quantum AI, "Quantum error correction below the surface code threshold", 2024 | 2408.13687 | 5 |
| Gidney, "Stim: a fast stabilizer circuit simulator", 2021 | 2103.02202 | 5 |
| Thaker et al., "Quantum Memory Hierarchies", 2006 | quant-ph/0604070 | 6 |
| Pino et al., "Demonstration of the trapped-ion quantum-CCD computer architecture", 2020 | 2003.01293 | 6 |
| Bluvstein et al., "A quantum processor based on coherent transport of entangled atom arrays", 2021 | 2112.03923 | 6 |
| Bluvstein et al., "Logical quantum processor based on reconfigurable atom arrays", 2023 | 2312.03982 | 6 |
| Jones et al., "Layered architecture for quantum computing", 2010 | 1010.5022 | 7 |
| Litinski, "A Game of Surface Codes", 2018 | 1808.02892 | 7 |
| Gidney, Ekerå, "How to factor 2048 bit RSA integers in 8 hours using 20 million noisy qubits", 2019 | 1905.09749 | 7 |
| Beverland et al., "Assessing requirements to scale to practical quantum advantage", 2022 | 2211.07629 | 7 |

Cited from memory, not yet verified:

- Knill, "Quantum computing with realistically noisy devices", Nature, 2005. Origin of the Pauli frame.
- Riesebos et al., "Pauli frames for quantum computer architectures", Design Automation Conference, 2017.
- Kielpinski, Monroe, Wineland, "Architecture for a large-scale ion-trap quantum computer", Nature, 2002.
- Van Meter and Horsman, "A blueprint for building a quantum computer", Communications of the ACM, 2013.

## Textbooks

- Nielsen and Chuang, *Quantum Computation and Quantum Information*. Chapter 8 covers Kraus operators and noise channels; chapter 10 covers error correction.
- Metodi, Faruque and Chong, *Quantum Computing for Computer Architects*. Closest textbook to the project's subject.
- Ding and Chong, *Quantum Computer Systems*. Compilation, qubit mapping and scheduling on current machines.
- De Wolf, *Quantum Computing: Lecture Notes*, arXiv 1907.09415. Grover and query complexity, rigorously.
