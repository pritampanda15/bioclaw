"use client";

import { JobList } from "@/components/jobs/JobList";

export default function JobsPage() {
  return (
    <div className="flex h-full flex-col">
      <div className="flex-1 overflow-y-auto p-6">
        <JobList />
      </div>
    </div>
  );
}
