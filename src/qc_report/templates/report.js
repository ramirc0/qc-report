(() => {
  const data = JSON.parse(document.getElementById("qc-data").textContent);
  const config = {
    displaylogo: false,
    responsive: true,
    modeBarButtonsToRemove: ["select2d", "lasso2d", "autoScale2d"],
  };

  function plot(host, figure, options = {}) {
    const layout = Object.assign({}, figure.layout, { template: data.template, autosize: true });
    Plotly.newPlot(host, figure.data, layout, Object.assign({}, config, options));
  }

  // Cards are small, so they skip the toolbar. Zoom and hover still work.
  document.querySelectorAll(".card .plot").forEach((host) =>
    plot(host, data.cards[host.dataset.card], { displayModeBar: false }),
  );

  // Sorting: numbers sort on data-value, sample names on their text. Blanks sort last.
  const table = document.querySelector(".samples-table");
  const body = table.tBodies[0];
  let sortedBy = null;
  let ascending = true;
  table.querySelectorAll("tr.cols th").forEach((th, column) => {
    th.addEventListener("click", () => {
      const text = th.dataset.type === "text";
      ascending = sortedBy === column ? !ascending : text;
      sortedBy = column;
      const key = (row) => {
        const cell = row.cells[column];
        if (text) return cell.textContent;
        const value = parseFloat(cell.dataset.value);
        return Number.isNaN(value) ? null : value;
      };
      const rows = Array.from(body.rows);
      rows.sort((a, b) => {
        const x = key(a);
        const y = key(b);
        if (x === null) return 1;
        if (y === null) return -1;
        return (x < y ? -1 : x > y ? 1 : 0) * (ascending ? 1 : -1);
      });
      body.append(...rows);
      table.querySelectorAll("th[aria-sort]").forEach((h) => h.removeAttribute("aria-sort"));
      th.setAttribute("aria-sort", ascending ? "ascending" : "descending");
    });
  });

  // Filtering: text, minimum level, outliers, and one audit flag.
  const search = document.getElementById("search");
  const shown = document.getElementById("shown");
  const flagFilter = document.getElementById("flag-filter");
  const filter = { minLevel: 0, outliers: false, flag: null };
  function apply() {
    const query = search.value.trim().toLowerCase();
    let count = 0;
    for (const row of body.rows) {
      const keep =
        row.dataset.sample.toLowerCase().includes(query) &&
        Number(row.dataset.level) >= filter.minLevel &&
        (!filter.outliers || Number(row.dataset.outliers) > 0) &&
        (!filter.flag || row.dataset.flags.split(" ").includes(filter.flag));
      row.hidden = !keep;
      count += keep;
    }
    shown.textContent = count;
  }
  search.addEventListener("input", apply);
  document.querySelectorAll(".filters button").forEach((button) => {
    button.addEventListener("click", () => {
      document.querySelectorAll(".filters button").forEach((b) => b.classList.remove("active"));
      button.classList.add("active");
      filter.outliers = button.hasAttribute("data-outliers");
      filter.minLevel = filter.outliers ? 0 : Number(button.dataset.minLevel);
      apply();
    });
  });
  function setFlag(row) {
    document.querySelectorAll(".flag-table tr").forEach((r) => r.classList.remove("active"));
    filter.flag = row && filter.flag !== row.dataset.flag ? row.dataset.flag : null;
    if (filter.flag) row.classList.add("active");
    flagFilter.hidden = !filter.flag;
    flagFilter.querySelector(".label").textContent = filter.flag ? row.dataset.label : "";
    apply();
  }
  document.querySelectorAll(".flag-table tr").forEach((row) => {
    row.addEventListener("click", () => {
      setFlag(row);
      document.getElementById("samples").scrollIntoView({ behavior: "smooth" });
    });
  });
  flagFilter.querySelector("button").addEventListener("click", () => setFlag(null));

  // Sample details open on #s-<sample>. Each loads samples/<sample>.js once.
  const requested = new Set();
  window.qcReport = {
    addSample({ sample, panels }) {
      const host = document.querySelector(`.panels[data-sample="${CSS.escape(sample)}"]`);
      host.replaceChildren();
      for (const panel of panels) {
        const figure = document.createElement("figure");
        figure.className = "panel";
        const title = document.createElement("h3");
        title.textContent = panel.title;
        figure.append(title);
        const content = document.createElement(panel.figure ? "div" : "img");
        if (panel.figure) {
          content.className = "plot";
          figure.append(content);
          host.append(figure);
          plot(content, panel.figure);
        } else {
          content.src = panel.image;
          content.alt = panel.title;
          figure.append(content);
          host.append(figure);
        }
      }
    },
  };
  function open() {
    const id = decodeURIComponent(location.hash.slice(1));
    if (!id.startsWith("s-")) return;
    const sample = id.slice(2);
    if (requested.has(sample)) return;
    requested.add(sample);
    const script = document.createElement("script");
    script.src = `${data.samplesDir}/${encodeURIComponent(sample)}.js`;
    script.onerror = () => {
      const host = document.querySelector(`.panels[data-sample="${CSS.escape(sample)}"]`);
      host.innerHTML = '<p class="empty">No figures for this sample.</p>';
    };
    document.body.append(script);
  }
  window.addEventListener("hashchange", open);
  document.addEventListener("keydown", (event) => {
    if (event.key === "Escape" && location.hash.startsWith("#s-")) location.hash = "samples";
  });
  open();
})();
