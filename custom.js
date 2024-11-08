document.addEventListener("DOMContentLoaded", function () {
  const headers = document.querySelectorAll(".book-toc li a");
  headers.forEach(header => {
    header.addEventListener("click", function () {
      const parentLi = header.parentElement;
      if (parentLi.classList.contains("active")) {
        parentLi.classList.remove("active");
      } else {
        parentLi.classList.add("active");
      }
    });
  });
});
